"""
    fit_jaiio_surrogates.jl

Fits `modelo_6b_refit` (Aproximado-2) and `modelo_7` (Aproximado-3) directly
against the new Julia `rigorous` V2 dataset (`JAIIO/data/v2_summary.csv`),
following the C_mix = 0.8 * C_tevap + 0.2 * C_bulk formulation from
surrogate_fit_common.py with Levenberg-Marquardt optimization + ForwardDiff
and 5-fold cross-validation.

Outputs:
  - JAIIO/results/fitted_surrogate_coeffs.csv
  - JAIIO/results/surrogate_fit_report.csv
"""

using CSV
using DataFrames
using ForwardDiff
using LinearAlgebra
using Printf
using QuadGK
using Random
using Statistics

include(joinpath(@__DIR__, "jaiio_common.jl"))
using .JAIIOCommon

const REGIME_CI = (
    (name=:low,  ci=0.495, pred = re -> re < 10.0),
    (name=:mid,  ci=0.496, pred = re -> 10.0 <= re < 100.0),
    (name=:high, ci=0.493, pred = re -> re >= 100.0),
)

# 8-point Gauss-Legendre nodes & weights on [0, 1] for exact ODE lifetime integration
# t_evap = (sqrt(pi) * d0^2 / (4 * dT)) * \int_0^1 (1 / C(sqrt(s)*d0)) ds
const _GL_NODES_RAW, _GL_WEIGHTS_RAW = QuadGK.gauss(8)
const GL_S = 0.5 .* (_GL_NODES_RAW .+ 1.0)
const GL_W = 0.5 .* _GL_WEIGHTS_RAW
const GL_SQRT_S = sqrt.(GL_S)

"""
Bounded Levenberg-Marquardt solver using exact ForwardDiff Jacobians.
Minimizes sum(abs2, residual_fn(theta)) subject to lb <= theta <= ub.
"""
function levenberg_marquardt(residual_fn, theta0::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64};
                             maxiter::Int=200, tol::Float64=1e-10)
    theta = clamp.(copy(theta0), lb, ub)
    r = residual_fn(theta)
    cost = dot(r, r)
    lambda = 1e-3

    for _ in 1:maxiter
        J = ForwardDiff.jacobian(residual_fn, theta)
        g = J' * r
        H = J' * J
        n = length(theta)

        accepted = false
        for _ in 1:15
            H_damped = H + lambda * Diagonal(max.(diag(H), 1e-8))
            step = -(H_damped \ g)
            theta_try = clamp.(theta + step, lb, ub)
            r_try = residual_fn(theta_try)
            cost_try = dot(r_try, r_try)
            if cost_try < cost
                rel_impr = (cost - cost_try) / max(cost, 1e-18)
                theta = theta_try
                r = r_try
                cost = cost_try
                lambda = max(lambda * 0.3, 1e-10)
                accepted = true
                if rel_impr < tol || norm(step) < tol
                    return theta, cost
                end
                break
            else
                lambda = min(lambda * 5.0, 1e10)
            end
        end
        if !accepted
            break
        end
    end
    return theta, cost
end

function prepare_fit_dataframe(v2_path::String)
    isfile(v2_path) || error("V2 summary CSV not found at $v2_path. Run run_jaiio_sweeps.jl first.")
    df = CSV.read(v2_path, DataFrame)
    rig = filter(r -> r.label == "rigorous" && r.t_evap_s < 29990.0 && isfinite(r.t_evap_s), df)
    rig.d0_mm = rig.d0_m .* 1000.0
    rig.dT = rig.T_inf_C .- rig.T_wb_C
    rig = filter(r -> r.dT > 1e-3, rig)
    rig.Theta = (rig.T0_C .- rig.T_wb_C) ./ rig.dT
    rig.C_bulk = rig.dV_dt_bulk_mm3s ./ (sqrt(π) .* rig.d0_mm .* rig.dT)
    rig.C_tevap = sqrt(π) .* (rig.d0_mm .^ 2) ./ (4.0 .* rig.dT .* rig.t_evap_s)
    rig.C_mix = 0.8 .* rig.C_tevap .+ 0.2 .* rig.C_bulk
    return rig
end

function model6_tevap_pred(p, d_mm, U, Tinf, dT, ci)
    A0, n0, A1, n1 = p
    t_ratio = max.(Tinf ./ 25.0, 1e-6)
    base_diff = A0 .* (t_ratio .^ n0)
    du_full   = (d_mm .* abs.(U)) .^ ci
    conv_pref = A1 .* (t_ratio .^ n1) .* du_full

    # Integrate 1 / C(sqrt(s)*d0) over s in [0, 1] via Gauss-Legendre quadrature
    invC_int = zeros(eltype(p), length(d_mm))
    for q in eachindex(GL_W)
        uq = GL_SQRT_S[q] ^ ci
        Cq = base_diff .+ conv_pref .* uq
        invC_int .+= GL_W[q] ./ Cq
    end
    return (sqrt(π) .* (d_mm .^ 2) ./ (4.0 .* dT)) .* invC_int
end

function model7_tevap_pred(p, d_mm, U, Tinf, HR, dT, Theta, ci)
    A0, n0, alpha, A1, n1, beta, gamma, delta = p
    t_ratio = max.(Tinf ./ 25.0, 1e-6)
    hr_term = clamp.(1.0 .- HR, 1e-4, 1.0)
    f_hr = hr_term .^ alpha
    g_hr = hr_term .^ beta
    h_th = @. exp(gamma * Theta + delta * (Theta ^ 2))

    base_diff = A0 .* (t_ratio .^ n0) .* f_hr .* h_th
    du_full   = (d_mm .* abs.(U)) .^ ci
    conv_pref = A1 .* (t_ratio .^ n1) .* du_full .* g_hr .* h_th

    invC_int = zeros(eltype(p), length(d_mm))
    for q in eachindex(GL_W)
        uq = GL_SQRT_S[q] ^ ci
        Cq = base_diff .+ conv_pref .* uq
        invC_int .+= GL_W[q] ./ Cq
    end
    return (sqrt(π) .* (d_mm .^ 2) ./ (4.0 .* dT)) .* invC_int
end

function fit_single_regime_m6b(sub::DataFrame, ci::Float64)
    d_mm   = Float64.(sub.d0_mm)
    U      = Float64.(sub.U_ms)
    Tinf   = Float64.(sub.T_inf_C)
    dT     = Float64.(sub.dT)
    t_true = Float64.(sub.t_evap_s)

    p0 = [5.5e-5, 0.24, 8.2e-5, 0.14]
    lb = [1e-6, -1.0, 1e-6, -1.0]
    ub = [5e-4,  2.0, 5e-4,  2.0]
    res_fn(p) = (model6_tevap_pred(p, d_mm, U, Tinf, dT, ci) .- t_true) ./ t_true
    popt, _ = levenberg_marquardt(res_fn, p0, lb, ub)
    return (popt[1], popt[2], popt[3], popt[4], ci)
end

function fit_single_regime_m7(sub::DataFrame, m6b_seed, ci::Float64)
    d_mm   = Float64.(sub.d0_mm)
    U      = Float64.(sub.U_ms)
    Tinf   = Float64.(sub.T_inf_C)
    HR     = Float64.(sub.HR)
    dT     = Float64.(sub.dT)
    Theta  = Float64.(sub.Theta)
    t_true = Float64.(sub.t_evap_s)

    p0 = [m6b_seed[1], m6b_seed[2], -0.04, m6b_seed[3], m6b_seed[4], -0.02, 0.005, 0.0]
    lb = [1e-6, -1.0, -1.5, 1e-6, -1.0, -1.5, -1.0, -1.0]
    ub = [5e-4,  2.0,  1.5, 5e-4,  2.0,  1.5,  1.0,  1.0]
    res_fn(p) = (model7_tevap_pred(p, d_mm, U, Tinf, HR, dT, Theta, ci) .- t_true) ./ t_true
    popt, _ = levenberg_marquardt(res_fn, p0, lb, ub)
    return (popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], ci)
end

function kfold_eval(sub::DataFrame, ci::Float64, is_m7::Bool; k::Int=5)
    n = nrow(sub)
    if n < k * 2
        return NaN
    end
    rng = MersenneTwister(42)
    perm = randperm(rng, n)
    fold_errs = Float64[]

    for fold in 1:k
        val_idx = perm[fold:k:n]
        trn_idx = setdiff(perm, val_idx)
        trn = sub[trn_idx, :]
        val = sub[val_idx, :]
        m6_p = fit_single_regime_m6b(trn, ci)
        yhat = if !is_m7
            model6_tevap_pred(collect(m6_p[1:4]), Float64.(val.d0_mm), Float64.(val.U_ms),
                              Float64.(val.T_inf_C), Float64.(val.dT), ci)
        else
            m7_p = fit_single_regime_m7(trn, m6_p, ci)
            model7_tevap_pred(collect(m7_p[1:8]), Float64.(val.d0_mm), Float64.(val.U_ms),
                              Float64.(val.T_inf_C), Float64.(val.HR), Float64.(val.dT),
                              Float64.(val.Theta), ci)
        end
        append!(fold_errs, abs.(yhat .- val.t_evap_s) ./ val.t_evap_s .* 100.0)
    end
    return mean(fold_errs)
end

function main()
    println("="^80)
    println("JAIIO 2026 SURROGATE RE-CALIBRATION AGAINST JULIA RIGOROUS PHYSICS")
    println("="^80)

    rig = prepare_fit_dataframe(joinpath(DATA_DIR, "v2_summary.csv"))
    @printf("Loaded %d rigorous reference cases from V2 dataset.\n", nrow(rig))

    m6_dict = Dict{Symbol, NTuple{5, Float64}}()
    m7_dict = Dict{Symbol, NTuple{9, Float64}}()
    report_rows = NamedTuple[]

    for reg in REGIME_CI
        sub = filter(r -> reg.pred(r.Re_0), rig)
        @printf("\nRegime %-4s (ci=%.3f): n = %d cases\n", string(reg.name), reg.ci, nrow(sub))

        p6 = fit_single_regime_m6b(sub, reg.ci)
        p7 = fit_single_regime_m7(sub, p6, reg.ci)
        m6_dict[reg.name] = p6
        m7_dict[reg.name] = p7

        # Full-sample & CV error on integrated t_evap_s
        y = Float64.(sub.t_evap_s)
        yhat6 = model6_tevap_pred(collect(p6[1:4]), Float64.(sub.d0_mm), Float64.(sub.U_ms),
                                  Float64.(sub.T_inf_C), Float64.(sub.dT), reg.ci)
        yhat7 = model7_tevap_pred(collect(p7[1:8]), Float64.(sub.d0_mm), Float64.(sub.U_ms),
                                  Float64.(sub.T_inf_C), Float64.(sub.HR), Float64.(sub.dT),
                                  Float64.(sub.Theta), reg.ci)

        mape6_train = mean(abs.(yhat6 .- y) ./ y .* 100.0)
        mape7_train = mean(abs.(yhat7 .- y) ./ y .* 100.0)
        mape6_cv    = kfold_eval(sub, reg.ci, false)
        mape7_cv    = kfold_eval(sub, reg.ci, true)

        @printf("  modelo_6b_refit: A0=%.6e, n0=%.5f, A1=%.6e, n1=%.5f | t_evap MAPE=%.2f%% (5-fold CV=%.2f%%)\n",
                p6[1], p6[2], p6[3], p6[4], mape6_train, mape6_cv)
        @printf("  modelo_7:        A0=%.6e, n0=%.5f, a=%.4f, A1=%.6e, n1=%.5f, b=%.4f, g=%.4f, d=%.4f | t_evap MAPE=%.2f%% (5-fold CV=%.2f%%)\n",
                p7[1], p7[2], p7[3], p7[4], p7[5], p7[6], p7[7], p7[8], mape7_train, mape7_cv)

        push!(report_rows, (model="modelo_6b_refit", regime=string(reg.name), n=nrow(sub), ci=reg.ci,
                            train_tevap_MAPE_pct=mape6_train, cv5_tevap_MAPE_pct=mape6_cv))
        push!(report_rows, (model="modelo_7", regime=string(reg.name), n=nrow(sub), ci=reg.ci,
                            train_tevap_MAPE_pct=mape7_train, cv5_tevap_MAPE_pct=mape7_cv))
    end

    m6b_coeffs = (low=m6_dict[:low], mid=m6_dict[:mid], high=m6_dict[:high])
    m7_coeffs  = (low=m7_dict[:low], mid=m7_dict[:mid], high=m7_dict[:high])

    coeffs_path = save_surrogate_coeffs(m6b_coeffs, m7_coeffs)
    report_path = joinpath(RESULTS_DIR, "surrogate_fit_report.csv")
    CSV.write(report_path, DataFrame(report_rows))

    println("\nSaved fitted coefficients to: ", coeffs_path)
    println("Saved fit CV report to:       ", report_path)
end

main()
