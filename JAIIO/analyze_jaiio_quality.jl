"""
    analyze_jaiio_quality.jl

Comprehensive Statistical, Physical-Regime, 2D ASABE, and Type-II ANCOVA
analysis suite for JAIIO 2026 & validation studies.

Reads:
  - JAIIO/data/v2_summary.csv (2,520 cases × 12 models)
  - JAIIO/data/v1_summary.csv (640 cases × 12 models, optional for 2D ASABE grids)

Writes tables to JAIIO/results/tables/:
  - paired_cases.csv
  - baseline_cases_with_t0_regimes.csv
  - global_metrics_all_outputs.csv
  - ranking_t_evap_s.csv
  - diffusion_vs_convection.csv
  - diffusion_by_d0.csv, diffusion_by_HR.csv, diffusion_by_Tinf.csv
  - convection_by_U.csv
  - t_evap_s_by_<regime>.csv (for Re_bin, d_bin, RH_bin, size_regime, humidity_regime,
                              re_regime_t0, sh_regime_t0, nu_regime_t0,
                              bm_regime_t0, bt_regime_t0, transport_regime_t0,
                              thermal_regime_t0, combined_regime_t0)
  - worst_cases_t_evap_s.csv (Top 25 worst cases per model)
  - top10_critical_by_model.csv (Top 10 critical errors per model for REPORT_JAIIO)
  - asabe_2d_v1_tables.csv (V1 PureDiff & Conv 2D ASABE × ThermalForce grids)
  - ancova_rigorous_tevap.csv, ancova_rigorous_dvdt.csv, ancova_error_<model>.csv
"""

using CSV
using DataFrames
using LinearAlgebra
using Printf
using Statistics

include(joinpath(@__DIR__, "jaiio_common.jl"))
using .JAIIOCommon

const INPUT_KEYS = [:d0_m, :U_ms, :T0_C, :T_inf_C, :HR]
const CANDIDATE_OUTPUTS = [
    :t_evap_s, :m_dot_0_kgs, :m_dot_median_kgs, :dV_dt_bulk_mm3s,
    :T_s0_C, :B_M_0, :B_T_0, :Re_0, :Nu_star_0, :Sh_star_0,
]

# -----------------------------------------------------------------------------
# Pure-Julia Statistical Metrics (Parity with model_quality_analysis_v6.py)
# -----------------------------------------------------------------------------

function quantile_linear(x::AbstractVector{Float64}, q::Float64)
    isempty(x) && return NaN
    return quantile(x, q)
end

function compute_paired_metrics(yt_raw::AbstractVector, yp_raw::AbstractVector)
    mask = isfinite.(yt_raw) .& isfinite.(yp_raw)
    yt = Float64.(yt_raw[mask])
    yp = Float64.(yp_raw[mask])
    n = length(yt)
    if n == 0
        return nothing
    end

    err = yp .- yt
    valid_rel = abs.(yt) .> 1e-15
    rel = fill(NaN, n)
    rel[valid_rel] = (err[valid_rel] ./ yt[valid_rel]) .* 100.0
    are = abs.(rel)
    rel_fin = filter(isfinite, rel)
    are_fin = filter(isfinite, are)

    mae = mean(abs.(err))
    rmse = sqrt(mean(abs2.(err)))
    span = maximum(yt) - minimum(yt)
    nrmse = span > 0 ? (rmse / span * 100.0) : (mean(abs.(yt)) > 0 ? rmse / mean(abs.(yt)) * 100.0 : NaN)

    mape = isempty(are_fin) ? NaN : mean(are_fin)
    medape = isempty(are_fin) ? NaN : median(are_fin)
    smape_denom = abs.(yt) .+ abs.(yp)
    sm_mask = smape_denom .> 1e-15
    smape = any(sm_mask) ? mean(200.0 .* abs.(err[sm_mask]) ./ smape_denom[sm_mask]) : NaN

    bias = mean(err)
    bias_pct = isempty(rel_fin) ? NaN : mean(rel_fin)
    medae = median(abs.(err))
    p95 = isempty(are_fin) ? NaN : quantile_linear(are_fin, 0.95)
    p99 = isempty(are_fin) ? NaN : quantile_linear(are_fin, 0.99)
    max_are = isempty(are_fin) ? NaN : maximum(are_fin)

    ybar = mean(yt)
    ss_res = sum(abs2.(err))
    ss_tot = sum(abs2.(yt .- ybar))
    r2 = ss_tot > 0 ? (1.0 - ss_res / ss_tot) : NaN

    ccc = if n >= 2
        mx, my = mean(yt), mean(yp)
        vx, vy = var(yt; corrected=false), var(yp; corrected=false)
        cov_xy = mean((yt .- mx) .* (yp .- my))
        denom = vx + vy + (mx - my)^2
        denom > 0 ? (2.0 * cov_xy / denom) : NaN
    else
        NaN
    end

    return (
        n=n, MAE=mae, RMSE=rmse, NRMSE_pct=nrmse,
        MAPE_pct=mape, MedAPE_pct=medape, sMAPE_pct=smape,
        Bias=bias, Bias_pct=bias_pct, MedAE=medae,
        P95_ARE_pct=p95, P99_ARE_pct=p99, Max_ARE_pct=max_are,
        R2=r2, CCC=ccc
    )
end

# -----------------------------------------------------------------------------
# Pure-Julia Incomplete Beta & F-distribution p-value for Type-II ANCOVA
# -----------------------------------------------------------------------------

# Continued fraction for regularized incomplete beta I_x(a, b) (Lentz method)
function _betacf(a::Float64, b::Float64, x::Float64)
    MAXIT = 200
    EPS = 3e-14
    FPMIN = 1e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    abs(d) < FPMIN && (d = FPMIN)
    d = 1.0 / d
    h = d
    for m in 1:MAXIT
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        abs(d) < FPMIN && (d = FPMIN)
        c = 1.0 + aa / c
        abs(c) < FPMIN && (c = FPMIN)
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        abs(d) < FPMIN && (d = FPMIN)
        c = 1.0 + aa / c
        abs(c) < FPMIN && (c = FPMIN)
        d = 1.0 / d
        del = d * c
        h *= del
        abs(del - 1.0) < EPS && break
    end
    return h
end

function _log_beta(a::Float64, b::Float64)
    # Stirling-accurate log Gamma via Lanczos
    return _lgamma(a) + _lgamma(b) - _lgamma(a + b)
end

function _lgamma(z::Float64)
    # Use Base / SpecialFunctions-free Lanczos g=7, n=9 approximation (accurate to ~1e-14)
    p = (
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    )
    if z < 0.5
        return log(π / sin(π * z)) - _lgamma(1.0 - z)
    end
    z -= 1.0
    x = p[1]
    for i in 2:9
        x += p[i] / (z + i - 1)
    end
    t = z + 7.5
    return 0.5 * log(2π) + (z + 0.5) * log(t) - t + log(x)
end

function regularized_beta(x::Float64, a::Float64, b::Float64)
    x <= 0.0 && return 0.0
    x >= 1.0 && return 1.0
    bt = exp(_lgamma(a + b) - _lgamma(a) - _lgamma(b) + a * log(x) + b * log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0)
        return bt * _betacf(a, b, x) / a
    else
        return 1.0 - bt * _betacf(b, a, 1.0 - x) / b
    end
end

function f_sf(F_stat::Float64, df1::Int, df2::Int)
    (!isfinite(F_stat) || F_stat <= 0.0 || df1 <= 0 || df2 <= 0) && return NaN
    x = df2 / (df2 + df1 * F_stat)
    return regularized_beta(x, df2 / 2.0, df1 / 2.0)
end

# -----------------------------------------------------------------------------
# Regime Definitions (Classic + t0-derived)
# -----------------------------------------------------------------------------

function sign_mag_quantile_labels(vals::AbstractVector{Float64})
    n = length(vals)
    out = fill("missing", n)
    neg_idx = findall(v -> isfinite(v) && v < -1e-8, vals)
    pos_idx = findall(v -> isfinite(v) && v >  1e-8, vals)
    nz_idx  = findall(v -> isfinite(v) && abs(v) <= 1e-8, vals)
    for i in nz_idx
        out[i] = "near_zero"
    end
    for (idxs, sgn) in ((neg_idx, "neg"), (pos_idx, "pos"))
        isempty(idxs) && continue
        mags = abs.(vals[idxs])
        q1, q2 = quantile(mags, [1/3, 2/3])
        for (loc, v) in zip(idxs, mags)
            out[loc] = v <= q1 ? "$(sgn)_low" : (v <= q2 ? "$(sgn)_mid" : "$(sgn)_high")
        end
    end
    return out
end

function build_baseline_t0_regimes(df::DataFrame)
    base = filter(r -> r.label == "rigorous", df)
    base = unique(base, INPUT_KEYS)
    base.d0_um = base.d0_m .* 1e6
    base.deltaT_inf_s0_C = base.T_inf_C .- base.T_s0_C
    base.abs_deltaT_inf_s0_C = abs.(base.deltaT_inf_s0_C)

    base.size_regime = [d < 100.0 ? "small_<100um" : (d < 500.0 ? "medium_100_500um" : "large_>=500um") for d in base.d0_um]
    base.humidity_regime = [h < 0.3 ? "dry_<0.3" : (h < 0.7 ? "intermediate_0.3_0.7" : "humid_>=0.7") for h in base.HR]
    base.re_regime_t0 = [re < 10.0 ? "Re<10" : (re < 100.0 ? "10<=Re<100" : "Re>=100") for re in base.Re_0]
    base.sh_regime_t0 = [s <= 2.5 ? "Sh*<=2.5" : (s <= 6.0 ? "2.5<Sh*<=6" : "Sh*>6") for s in base.Sh_star_0]
    base.nu_regime_t0 = [nu <= 2.5 ? "Nu*<=2.5" : (nu <= 6.0 ? "2.5<Nu*<=6" : "Nu*>6") for nu in base.Nu_star_0]
    base.bm_regime_t0 = sign_mag_quantile_labels(Float64.(base.B_M_0))
    base.bt_regime_t0 = sign_mag_quantile_labels(Float64.(base.B_T_0))

    base.transport_regime_t0 = map(base.Re_0, base.Sh_star_0) do re, sh
        if re < 10.0 && sh <= 2.5
            "diffusive"
        elseif (10.0 <= re < 100.0) || (2.5 < sh <= 6.0)
            "mixed"
        elseif re >= 100.0 || sh > 6.0
            "convective_enhanced"
        else
            "unclassified"
        end
    end

    dt_vals = filter(isfinite, base.abs_deltaT_inf_s0_C)
    dt_q = isempty(dt_vals) ? NaN : quantile(dt_vals, 2/3)
    base.thermal_regime_t0 = map(base.B_T_0, base.abs_deltaT_inf_s0_C) do bt, dt
        if bt <= 0.0
            "cooling_limited"
        elseif dt <= dt_q
            "thermally_balanced"
        else
            "strong_thermal_gradient"
        end
    end

    base.combined_regime_t0 = base.size_regime .* " | " .* base.humidity_regime .* " | " .*
                              base.transport_regime_t0 .* " | " .* base.thermal_regime_t0
    return base
end

function build_paired_dataset(df::DataFrame)
    base_regimes = build_baseline_t0_regimes(df)
    ref = filter(r -> r.label == "rigorous", df)
    ref = unique(ref, INPUT_KEYS)

    # Rename reference outputs to <col>__ref
    ref_keep = unique(vcat(INPUT_KEYS, [:Re_0], CANDIDATE_OUTPUTS))
    ref_sub = select(ref, [c for c in ref_keep if hasproperty(ref, c)])
    for c in vcat(CANDIDATE_OUTPUTS, [:Re_0])
        if hasproperty(ref_sub, c) && !hasproperty(ref_sub, Symbol("$(c)__ref"))
            rename!(ref_sub, c => Symbol("$(c)__ref"))
        end
    end

    # Join t0 regime columns onto ref_sub
    t0_cols = [:size_regime, :humidity_regime, :re_regime_t0, :sh_regime_t0, :nu_regime_t0,
               :bm_regime_t0, :bt_regime_t0, :transport_regime_t0, :thermal_regime_t0,
               :combined_regime_t0, :deltaT_inf_s0_C, :abs_deltaT_inf_s0_C]
    reg_join = select(base_regimes, vcat(INPUT_KEYS, t0_cols))
    ref_full = innerjoin(ref_sub, reg_join; on=INPUT_KEYS)

    # Pair every non-rigorous model with ref_full
    models = [m.id for m in MODEL_SPECS if m.id !== :rigorous && string(m.id) in df.label]
    paired_parts = DataFrame[]
    for m_id in models
        m_str = string(m_id)
        mod_df = filter(r -> r.label == m_str, df)
        m_keep = vcat(INPUT_KEYS, [c for c in CANDIDATE_OUTPUTS if hasproperty(mod_df, c)])
        mod_sub = select(mod_df, m_keep)
        merged = innerjoin(mod_sub, ref_full; on=INPUT_KEYS)
        merged.model .= m_str
        merged.display_name .= display_name(m_id)
        push!(paired_parts, merged)
    end

    paired = vcat(paired_parts...)
    paired.d0_um = paired.d0_m .* 1e6
    paired.Re_for_regime = paired.Re_0__ref

    # Classic bins (note cut-point semantics match pd.cut include_lowest=True)
    paired.Re_bin = [re <= 10.0 ? "Re<=10" : (re <= 100.0 ? "10<Re<=100" : "Re>100") for re in paired.Re_for_regime]
    paired.d_bin  = [d <= 100.0 ? "d<=100 um" : (d <= 500.0 ? "100<d<=500 um" : "d>500 um") for d in paired.d0_um]
    paired.RH_bin = [h <= 0.3 ? "RH<=0.3" : (h <= 0.7 ? "0.3<RH<=0.7" : "RH>0.7") for h in paired.HR]

    # Signed relative error and absolute relative error on t_evap_s
    paired.rel_err_pct = ((paired.t_evap_s .- paired.t_evap_s__ref) ./ paired.t_evap_s__ref) .* 100.0
    paired.are_pct = abs.(paired.rel_err_pct)

    return paired, base_regimes
end

# -----------------------------------------------------------------------------
# Type-II ANCOVA Engine (Pure Julia QR-OLS)
# -----------------------------------------------------------------------------

function run_ancova_type2(df_in::DataFrame, response_col::Symbol; tag::String)
    df = copy(df_in)
    df = filter(r -> isfinite(r[response_col]), df)
    y = Float64.(df[!, response_col])
    n = length(y)

    # Categoricals: d0_um_cat (drop ref level), HR_cat (drop ref level)
    d0_cats = sort(unique(round.(df.d0_um; digits=1)))
    hr_cats = sort(unique(round.(Int, df.HR .* 100)))

    # Build term specs: (name, matrix_columns)
    terms = Tuple{String, Matrix{Float64}}[]
    if length(d0_cats) > 1
        M_d0 = zeros(n, length(d0_cats) - 1)
        for j in 2:length(d0_cats)
            M_d0[:, j-1] = Float64.(round.(df.d0_um; digits=1) .== d0_cats[j])
        end
        push!(terms, ("C(d0_um)", M_d0))
    end
    if length(hr_cats) > 1
        M_hr = zeros(n, length(hr_cats) - 1)
        for j in 2:length(hr_cats)
            M_hr[:, j-1] = Float64.(round.(Int, df.HR .* 100) .== hr_cats[j])
        end
        push!(terms, ("C(HR)", M_hr))
    end
    for cov in (:U_ms, :T_inf_C, :T0_C, :Theta, :B_M_0)
        if hasproperty(df, cov) && all(isfinite, df[!, cov])
            v = Float64.(df[!, cov])
            if std(v) > 1e-12
                push!(terms, (string(cov), reshape(v, n, 1)))
            end
        end
    end

    # Assemble full design matrix [1, terms...]
    X_full = hcat(ones(n, 1), [t[2] for t in terms]...)
    beta_full = X_full \ y
    resid_full = y .- X_full * beta_full
    sse_full = dot(resid_full, resid_full)
    p_full = size(X_full, 2)
    df_resid = n - p_full
    mse_resid = sse_full / max(df_resid, 1)

    rows = NamedTuple[]
    # Compute Type-II SS by dropping each term one at a time from full model
    col_offset = 2
    sum_ss_effects = 0.0
    for (tname, M_t) in terms
        df_t = size(M_t, 2)
        keep_cols = vcat(1:(col_offset - 1), (col_offset + df_t):p_full)
        X_red = X_full[:, keep_cols]
        beta_red = X_red \ y
        sse_red = sum(abs2, y .- X_red * beta_red)
        ss_t = max(sse_red - sse_full, 0.0)
        sum_ss_effects += ss_t
        col_offset += df_t

        F_val = (ss_t / df_t) / max(mse_resid, 1e-30)
        pval = f_sf(F_val, df_t, df_resid)
        push!(rows, (term=tname, sum_sq=ss_t, df=df_t, F=F_val, p_value=pval, eta_sq=NaN))
    end
    push!(rows, (term="Residual", sum_sq=sse_full, df=df_resid, F=NaN, p_value=NaN, eta_sq=NaN))

    ss_total_type2 = sum_ss_effects + sse_full
    out_df = DataFrame(rows)
    out_df.eta_sq = out_df.sum_sq ./ max(ss_total_type2, 1e-30)

    out_path = joinpath(TABLES_DIR, "ancova_$(tag).csv")
    CSV.write(out_path, out_df)
    return out_df
end

# -----------------------------------------------------------------------------
# V1 2D ASABE × ThermalForce Grids
# -----------------------------------------------------------------------------

asabe_cat(d0_m) = begin
    d_um = d0_m * 1e6
    d_um < 60  ? "1_XF" :
    d_um < 145 ? "2_VF" :
    d_um < 225 ? "3_F"  :
    d_um < 325 ? "4_M"  :
    d_um < 400 ? "5_C"  :
    d_um < 500 ? "6_VC" :
    d_um < 650 ? "7_XC" : "8_UC"
end

thermal_cat(dT) = dT < 3 ? "1_VeryLow" : (dT < 10 ? "2_Low" : (dT < 20 ? "3_Medium" : (dT < 30 ? "4_High" : "5_VeryHigh")))

function analyze_v1_2d_asabe(v1_path::String)
    if !isfile(v1_path)
        println("V1 dataset not found at $v1_path; skipping 2D ASABE analysis.")
        return DataFrame()
    end
    df = CSV.read(v1_path, DataFrame)
    rig = filter(r -> r.label == "rigorous" && isfinite(r.t_evap_s), df)
    rig = unique(rig, INPUT_KEYS)
    rig.ASABE = asabe_cat.(rig.d0_m)
    rig.ThermalForce = thermal_cat.(rig.T_inf_C .- rig.T_wb_C)
    rig_ref = select(rig, vcat(INPUT_KEYS, [:t_evap_s, :ASABE, :ThermalForce]))
    rename!(rig_ref, :t_evap_s => :t_evap_s__ref)

    rows = NamedTuple[]
    models = [m.id for m in MODEL_SPECS if m.id !== :rigorous && string(m.id) in df.label]
    for m_id in models
        m_str = string(m_id)
        mod_df = filter(r -> r.label == m_str, df)
        mrg = innerjoin(select(mod_df, vcat(INPUT_KEYS, [:t_evap_s])), rig_ref; on=INPUT_KEYS)
        mrg.are_pct = abs.(mrg.t_evap_s .- mrg.t_evap_s__ref) ./ mrg.t_evap_s__ref .* 100.0
        for (reg_name, sub) in (("PureDiff", filter(r -> r.U_ms == 0.0, mrg)),
                                ("Conv",     filter(r -> r.U_ms > 0.0, mrg)))
            nrow(sub) == 0 && continue
            g = combine(groupby(sub, [:ASABE, :ThermalForce]),
                :are_pct => mean => :mean_are_pct,
                :are_pct => (x -> length(x) > 1 ? std(x) : 0.0) => :std_are_pct,
                :are_pct => maximum => :max_are_pct,
                nrow => :n
            )
            g.model .= m_str
            g.display_name .= display_name(m_id)
            g.regime .= reg_name
            for r in eachrow(g)
                push!(rows, (; model=r.model, display_name=r.display_name, regime=r.regime,
                               ASABE=r.ASABE, ThermalForce=r.ThermalForce,
                               n=r.n, mean_are_pct=r.mean_are_pct,
                               std_are_pct=r.std_are_pct, max_are_pct=r.max_are_pct))
            end
        end
    end

    out_df = DataFrame(rows)
    CSV.write(joinpath(TABLES_DIR, "asabe_2d_v1_tables.csv"), out_df)
    return out_df
end

# -----------------------------------------------------------------------------
# Main Analysis Driver
# -----------------------------------------------------------------------------

function main()
    println("="^80)
    println("JAIIO 2026 QUALITY, REGIME, 2D ASABE & ANCOVA ANALYSIS (Julia)")
    println("="^80)

    v2_path = joinpath(DATA_DIR, "v2_summary.csv")
    isfile(v2_path) || error("Missing $v2_path. Run run_jaiio_sweeps.jl first.")
    df_v2 = CSV.read(v2_path, DataFrame)
    df_v2_clean = filter(r -> r.t_evap_s < 29990.0 && isfinite(r.t_evap_s), df_v2)

    paired, base_t0 = build_paired_dataset(df_v2_clean)
    CSV.write(joinpath(TABLES_DIR, "baseline_cases_with_t0_regimes.csv"), base_t0)
    CSV.write(joinpath(TABLES_DIR, "paired_cases.csv"), paired)
    @printf("Paired V2 dataset: %d rows (%d unique rigorous baseline conditions)\n", nrow(paired), nrow(base_t0))

    # 1. Global metrics across all candidate outputs
    order_idx = Dict(string(m.id) => i for (i, m) in enumerate(MODEL_SPECS))
    global_rows = NamedTuple[]
    for m_str in unique(paired.model)
        sdf = filter(r -> r.model == m_str, paired)
        for out_col in CANDIDATE_OUTPUTS
            ref_col = Symbol("$(out_col)__ref")
            (hasproperty(sdf, out_col) && hasproperty(sdf, ref_col)) || continue
            met = compute_paired_metrics(sdf[!, ref_col], sdf[!, out_col])
            met === nothing && continue
            push!(global_rows, (; model=m_str, display_name=display_name(m_str), output=string(out_col), met...))
        end
    end
    df_global = DataFrame(global_rows)
    CSV.write(joinpath(TABLES_DIR, "global_metrics_all_outputs.csv"), df_global)

    # 2. Primary ranking on t_evap_s
    rank_df = filter(r -> r.output == "t_evap_s", df_global)
    rank_df.abs_bias = abs.(rank_df.Bias)
    sort!(rank_df, [:MAPE_pct, :RMSE, :P95_ARE_pct, :abs_bias])
    rank_df.rank = 1:nrow(rank_df)
    select!(rank_df, Not(:abs_bias))
    CSV.write(joinpath(TABLES_DIR, "ranking_t_evap_s.csv"), rank_df)
    println("\n--- Global V2 Ranking (t_evap_s) ---")
    println(rank_df[:, [:rank, :model, :display_name, :n, :MAPE_pct, :MedAPE_pct, :Bias_pct, :P95_ARE_pct, :Max_ARE_pct, :R2, :CCC]])

    # 3. Regime tables (classic bins + t0 regimes)
    regime_cols = [
        :Re_bin, :d_bin, :RH_bin,
        :size_regime, :humidity_regime, :re_regime_t0, :sh_regime_t0, :nu_regime_t0,
        :bm_regime_t0, :bt_regime_t0, :transport_regime_t0, :thermal_regime_t0, :combined_regime_t0
    ]
    for rcol in regime_cols
        rows = NamedTuple[]
        for m_str in unique(paired.model)
            mdf = filter(r -> r.model == m_str, paired)
            for grp in groupby(mdf, rcol)
                met = compute_paired_metrics(grp.t_evap_s__ref, grp.t_evap_s)
                met === nothing && continue
                push!(rows, (; model=m_str, display_name=display_name(m_str),
                               regime_val=string(grp[1, rcol]), met...))
            end
        end
        reg_df = DataFrame(rows)
        rename!(reg_df, :regime_val => rcol)
        CSV.write(joinpath(TABLES_DIR, "t_evap_s_by_$(rcol).csv"), reg_df)
    end

    # 4. Pure Diffusion (U = 0) vs Forced Convection (U > 0) analysis
    dc_rows = NamedTuple[]
    for (reg_name, pred) in (("Pure Diffusion (U=0)", r -> r.U_ms == 0.0),
                             ("Forced Convection (U>0)", r -> r.U_ms > 0.0))
        sub_reg = filter(pred, paired)
        for m_str in unique(sub_reg.model)
            mdf = filter(r -> r.model == m_str, sub_reg)
            met = compute_paired_metrics(mdf.t_evap_s__ref, mdf.t_evap_s)
            met === nothing && continue
            push!(dc_rows, (; regime=reg_name, model=m_str, display_name=display_name(m_str), met...))
        end
    end
    df_dc = DataFrame(dc_rows)
    CSV.write(joinpath(TABLES_DIR, "diffusion_vs_convection.csv"), df_dc)

    # Detailed diffusion breakdowns (by d0_um, HR, T_inf_C) & convection breakdown (by U_ms)
    diff_df = filter(r -> r.U_ms == 0.0, paired)
    conv_df = filter(r -> r.U_ms > 0.0, paired)
    for (sub_df, grp_col, fname) in ((diff_df, :d0_um,   "diffusion_by_d0.csv"),
                                     (diff_df, :HR,      "diffusion_by_HR.csv"),
                                     (diff_df, :T_inf_C, "diffusion_by_Tinf.csv"),
                                     (conv_df, :U_ms,    "convection_by_U.csv"))
        rows = NamedTuple[]
        for m_str in unique(sub_df.model)
            mdf = filter(r -> r.model == m_str, sub_df)
            for grp in groupby(mdf, grp_col)
                ares = grp.are_pct
                push!(rows, (; model=m_str, display_name=display_name(m_str),
                               bin_val=grp[1, grp_col], n=nrow(grp),
                               MAPE_pct=mean(ares),
                               Std_ARE_pct=(nrow(grp) > 1 ? std(ares) : 0.0),
                               Max_ARE_pct=maximum(ares),
                               Bias_pct=mean(grp.rel_err_pct)))
            end
        end
        CSV.write(joinpath(TABLES_DIR, fname), DataFrame(rows))
    end

    # 5. Worst Cases (Top 25) and Critical Errors (Top 10 for REPORT_JAIIO)
    worst_parts = DataFrame[]
    top10_parts = DataFrame[]
    for m_str in unique(paired.model)
        mdf = filter(r -> r.model == m_str, paired)
        sort!(mdf, :are_pct, rev=true)
        keep_cols = [:model, :display_name, :d0_m, :d0_um, :U_ms, :T0_C, :T_inf_C, :HR,
                     :t_evap_s__ref, :t_evap_s, :rel_err_pct, :are_pct,
                     :Re_for_regime, :size_regime, :humidity_regime, :transport_regime_t0, :thermal_regime_t0]
        push!(worst_parts, first(select(mdf, keep_cols), min(25, nrow(mdf))))
        push!(top10_parts, first(select(mdf, keep_cols), min(10, nrow(mdf))))
    end
    CSV.write(joinpath(TABLES_DIR, "worst_cases_t_evap_s.csv"), vcat(worst_parts...))
    CSV.write(joinpath(TABLES_DIR, "top10_critical_by_model.csv"), vcat(top10_parts...))

    # 6. V1 2D ASABE Analysis
    analyze_v1_2d_asabe(joinpath(DATA_DIR, "v1_summary.csv"))

    # 7. Type-II ANCOVA on rigorous physics & surrogate relative errors
    rig_anc = filter(r -> r.label == "rigorous", df_v2_clean)
    rig_anc.d0_um = rig_anc.d0_m .* 1e6
    rig_anc.dT = rig_anc.T_inf_C .- rig_anc.T_wb_C
    rig_anc.Theta = [dt > 1e-6 ? (t0 - twb) / dt : 0.0 for (t0, twb, dt) in zip(rig_anc.T0_C, rig_anc.T_wb_C, rig_anc.dT)]

    run_ancova_type2(rig_anc, :t_evap_s; tag="rigorous_tevap")
    run_ancova_type2(rig_anc, :dV_dt_bulk_mm3s; tag="rigorous_dvdt")

    for m_anc in ("modelo_6b_refit", "modelo_7", "modelo_2b", "level_2", "uniform", "iso_twb")
        sub_err = filter(r -> r.model == m_anc, paired)
        nrow(sub_err) < 20 && continue
        sub_err = innerjoin(sub_err, select(rig_anc, vcat(INPUT_KEYS, [:Theta])); on=INPUT_KEYS)
        # Ensure ANCOVA uses rigorous baseline B_M_0 (surrogates have NaN for B_M_0)
        if hasproperty(sub_err, :B_M_0__ref)
            sub_err.B_M_0 = sub_err.B_M_0__ref
        end
        run_ancova_type2(sub_err, :rel_err_pct; tag="error_$(m_anc)")
    end

    println("\nAll analysis CSV tables written to: ", TABLES_DIR)
end

main()
