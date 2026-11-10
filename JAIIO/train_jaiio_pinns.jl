"""
    train_jaiio_pinns.jl

Generates/caches rigorous collocation trajectory datasets (V2: 2,520 cases; V1: 640 cases),
trains and cross-validates four neural & Physics-Informed architectures in native Julia
(Zygote + Optimisers.AdamW):
  1. PINN-1Eq     (:pinn_c)      - Physics-Informed residual d²-law + monotonicity penalties
  2. PINN-2Eq     (:pinn_2eq)    - Coupled [V(t), T_bar(t)] + Neural Biot Ts + AS energy balance
  3. PI-DeepONet  (:pi_deeponet) - Solver-free operator network + hard BCs + ODE collocation
  4. ResMLP-Puro  (:resmlp)      - Pure data-driven MLP baseline (ablation control)

Outputs:
  - JAIIO/data/pinn_v2_trajectories.npz
  - JAIIO/data/pinn_v1_trajectories.npz
  - JAIIO/results/pinn_weights.npz
  - JAIIO/results/pinn_training_history.csv
  - JAIIO/results/pinn_cv_metrics.csv
  - JAIIO/results/pinn_v1_ood_metrics.csv

Usage:
  julia --threads auto --project=julia JAIIO/train_jaiio_pinns.jl [--smoke] [--regen-data]
"""

using ArgParse
using CSV
using DataFrames
using Dates
using LinearAlgebra
using NPZ
using Optimisers
using OrdinaryDiffEq
using Printf
using QuadGK
using Random
using Statistics
using Zygote

include(joinpath(@__DIR__, "jaiio_common.jl"))
using .JAIIOCommon
using Evaporacion.Common
using Evaporacion.Droplet
using Evaporacion.Surrogates
using Evaporacion.Neural.PINNSurrogates

const TAU_GRID = Float64[
    0.0, 0.002, 0.005, 0.01, 0.02, 0.035, 0.05, 0.08, 0.12, 0.18,
    0.25, 0.33, 0.42, 0.50, 0.58, 0.66, 0.74, 0.81, 0.87, 0.92,
    0.95, 0.97, 0.985, 0.995
]
const N_TAU = length(TAU_GRID)

function parse_cli()
    s = ArgParseSettings(description="Train & cross-validate JAIIO PINN and Neural Operator models.")
    @add_arg_table! s begin
        "--smoke"
            help = "Fast smoke run (120 V2 cases, 60 V1 cases, short epochs)"
            action = :store_true
        "--regen-data"
            help = "Force regeneration of trajectory collocation caches"
            action = :store_true
    end
    return parse_args(s)
end

# =============================================================================
# 1. RIGOROUS COLLOCATION TRAJECTORY DATASET GENERATION
# =============================================================================

function sample_rigorous_trajectories(name::String, conditions)
    N = length(conditions)
    @printf("[%s] Generating rigorous (:brent) trajectories for %d cases (%d threads)...\n",
            name, N, Threads.nthreads())

    # Warmup JIT
    c1 = conditions[1]
    p1 = JAIIOCommon.build_physical_params(:rigorous, c1.d0_m, c1.U_ms, c1.T0_C, c1.T_inf_C, c1.HR; ts_solver=:brent)
    m01 = liquid_density_water(c1.T0_C + T_ZERO_C) * (π / 6.0) * (c1.d0_m^3)
    solve_droplet(p1, [m01, c1.T0_C], (0.0, 30000.0);
                  solver=AutoTsit5(Rosenbrock23(autodiff=false)),
                  atol=[min(1e-18, m01 * 1e-9), 1e-9], rtol=1e-6, maxiters=10^6)

    d0_m      = zeros(Float64, N)
    U_ms      = zeros(Float64, N)
    T0_C      = zeros(Float64, N)
    T_inf_C   = zeros(Float64, N)
    HR        = zeros(Float64, N)
    T_wb_C    = zeros(Float64, N)
    t_evap_s  = zeros(Float64, N)

    V_mm3  = zeros(Float64, N_TAU, N)
    d_mm   = zeros(Float64, N_TAU, N)
    T_bar  = zeros(Float64, N_TAU, N)
    Ts     = zeros(Float64, N_TAU, N)
    m_dot  = zeros(Float64, N_TAU, N)
    dV_dt  = zeros(Float64, N_TAU, N)
    dT_dt  = zeros(Float64, N_TAU, N)

    done_cnt = Threads.Atomic{Int}(0)
    t0 = time()

    Threads.@threads for i in 1:N
        c = conditions[i]
        d0_m[i]    = c.d0_m
        U_ms[i]    = c.U_ms
        T0_C[i]    = c.T0_C
        T_inf_C[i] = c.T_inf_C
        HR[i]      = c.HR
        twb        = stull_twb(c.T_inf_C, c.HR * 100.0)
        T_wb_C[i]  = twb

        params = JAIIOCommon.build_physical_params(:rigorous, c.d0_m, c.U_ms, c.T0_C, c.T_inf_C, c.HR; ts_solver=:brent)
        rl0 = liquid_density_water(c.T0_C + T_ZERO_C)
        m0  = rl0 * (π / 6.0) * (c.d0_m^3)
        sol = solve_droplet(params, [m0, c.T0_C], (0.0, 30000.0);
                            solver=AutoTsit5(Rosenbrock23(autodiff=false)),
                            atol=[min(1e-18, m0 * 1e-9), 1e-9], rtol=1e-6, maxiters=10^6)
        tev = Float64(sol.t[end])
        t_evap_s[i] = tev

        for k in 1:N_TAU
            tk = clamp(TAU_GRID[k] * tev, 0.0, tev)
            uk = sol(tk)
            mk = max(Float64(uk[1]), 1e-15)
            Tk = Float64(uk[2])
            Tsk = Float64(solve_surface_temperature(mk, Tk, params))
            alg = calc_droplet_algebraic(mk, Tk, Tsk, params)

            Vk = (mk / alg.rho_l) * 1e9
            dk = alg.d * 1e3
            dVdt_k = -(alg.m_dot / alg.rho_l) * 1e9

            term1 = 3.0 * alg.h * (c.T_inf_C - Tsk) / (alg.rho_l * alg.cp_l * alg.r_s)
            term2 = 3.0 * alg.r_dot_s * alg.DHev / (alg.r_s * alg.cp_l)
            term3 = 3.0 * alg.r_dot_s * (Tsk - Tk) / alg.r_s
            dTdt_k = term1 + term2 + term3

            V_mm3[k, i] = Vk
            d_mm[k, i]  = dk
            T_bar[k, i] = Tk
            Ts[k, i]    = Tsk
            m_dot[k, i] = Float64(alg.m_dot)
            dV_dt[k, i] = dVdt_k
            dT_dt[k, i] = dTdt_k
        end

        nd = Threads.atomic_add!(done_cnt, 1) + 1
        if nd % max(1, N ÷ 8) == 0 || nd == N
            @printf("  [%s] Trajectories: %d / %d (%.1f%%) [%.1fs]\n",
                    name, nd, N, 100.0 * nd / N, time() - t0)
        end
    end

    return Dict{String, Any}(
        "tau" => collect(TAU_GRID),
        "d0_m" => d0_m, "U_ms" => U_ms, "T0_C" => T0_C, "T_inf_C" => T_inf_C,
        "HR" => HR, "T_wb_C" => T_wb_C, "t_evap_s" => t_evap_s,
        "V_mm3" => V_mm3, "d_mm" => d_mm, "T_bar" => T_bar, "Ts" => Ts,
        "m_dot" => m_dot, "dV_dt" => dV_dt, "dT_dt" => dT_dt
    )
end

function build_or_load_trajectories(name::String, conditions, path::String; force_regen::Bool=false)
    if !force_regen && isfile(path)
        d = npzread(path)
        if length(d["t_evap_s"]) == length(conditions)
            @printf("[%s] Loaded cached trajectories (%d cases) from %s\n", name, length(conditions), path)
            return d
        end
    end
    d = sample_rigorous_trajectories(name, conditions)
    tmp = path * ".tmp.npz"
    npzwrite(tmp, d)
    mv(tmp, path; force=true)
    @printf("[%s] Saved trajectory cache to %s\n", name, path)
    return d
end

function slice_dataset(d::Dict{String, Any}, idx::AbstractVector{Int})
    out = Dict{String, Any}("tau" => d["tau"])
    for k in ("d0_m", "U_ms", "T0_C", "T_inf_C", "HR", "T_wb_C", "t_evap_s")
        out[k] = d[k][idx]
    end
    for k in ("V_mm3", "d_mm", "T_bar", "Ts", "m_dot", "dV_dt", "dT_dt")
        out[k] = d[k][:, idx]
    end
    return out
end

# Helper to compute standard metrics cleanly
@inline mape_pct(yhat, y) = mean(abs.(yhat .- y) ./ max.(abs.(y), 1e-18)) * 100.0
@inline mae_val(yhat, y)  = mean(abs.(yhat .- y))
@inline max_are_pct(yhat, y) = maximum(abs.(yhat .- y) ./ max.(abs.(y), 1e-18)) * 100.0

@inline function _l2_params(ps::NamedTuple)
    s = 0.0
    for l in ps.layers
        s += sum(abs2, l.W)
    end
    return s
end

# =============================================================================
# 2. TRAIN MODEL 1: PINN-C (1-Eq Physics-Informed Residual d²-Law)
# =============================================================================

function prepare_pinn_c_tensors(d::Dict{String, Any})
    N = length(d["t_evap_s"])
    nq = length(PINNSurrogates.GL_S) # 8
    M = N * nq

    X_gl      = zeros(Float64, 9, M)
    C_base_gl = zeros(Float64, M)
    X_u       = zeros(Float64, 9, M)
    C_base_u  = zeros(Float64, M)
    X_hr      = zeros(Float64, 9, M)
    C_base_hr = zeros(Float64, M)
    prefactor = zeros(Float64, N)
    t_true    = Float64.(d["t_evap_s"])

    for i in 1:N
        d0_mm = d["d0_m"][i] * 1e3
        U     = d["U_ms"][i]
        T0    = d["T0_C"][i]
        Tinf  = d["T_inf_C"][i]
        HR    = d["HR"][i]
        sc    = PINNSurrogates.compute_static_thermo_scales(Tinf, HR, T0)
        dT    = max(sc.dT, 1e-3)
        prefactor[i] = sqrt(π) * (d0_mm^2) / (4.0 * dT)

        # Perturbed conditions for physics monotonicity collocation
        U_p  = U + 0.15
        HR_p = min(HR + 0.05, 0.98)
        sc_hr = PINNSurrogates.compute_static_thermo_scales(Tinf, HR_p, T0)

        for q in 1:nq
            col = (i - 1) * nq + q
            dq  = d0_mm * PINNSurrogates.GL_SQRT_S[q]

            X_gl[:, col]  = collect(pinn_c_features(dq, d0_mm, U, T0, Tinf, HR, sc.Theta0, sc.BM_wb))
            C_base_gl[col] = calc_C_model7_smooth(dq, U, Tinf, HR, sc.T_wb, T0)

            X_u[:, col]   = collect(pinn_c_features(dq, d0_mm, U_p, T0, Tinf, HR, sc.Theta0, sc.BM_wb))
            C_base_u[col] = calc_C_model7_smooth(dq, U_p, Tinf, HR, sc.T_wb, T0)

            X_hr[:, col]  = collect(pinn_c_features(dq, d0_mm, U, T0, Tinf, HR_p, sc_hr.Theta0, sc_hr.BM_wb))
            C_base_hr[col] = calc_C_model7_smooth(dq, U, Tinf, HR_p, sc_hr.T_wb, T0)
        end
    end

    return (
        N=N, nq=nq, X_gl=X_gl, C_base_gl=C_base_gl,
        X_u=X_u, C_base_u=C_base_u, X_hr=X_hr, C_base_hr=C_base_hr,
        prefactor=prefactor, t_true=t_true,
        gl_w=PINNSurrogates.GL_W
    )
end

function pinn_c_eval_batch(ps, tens)
    out = vec(mlp_forward(ps, tens.X_gl))
    C_gl = tens.C_base_gl .* exp.(0.15 .* tanh.(out))
    C_mat = reshape(C_gl, tens.nq, tens.N)
    invC = vec(sum(tens.gl_w ./ C_mat; dims=1))
    t_hat = tens.prefactor .* invC
    return t_hat, C_gl
end

function train_pinn_c(train_d, val_d; epochs::Int=800, lr::Float64=3e-3, seed::Int=42, log_history::Bool=false)
    rng = MersenneTwister(seed)
    ps = init_mlp_params(rng, [9, 24, 24, 1]; out_scale=0.01)
    tr = prepare_pinn_c_tensors(train_d)
    va = val_d === nothing ? nothing : prepare_pinn_c_tensors(val_d)

    opt = Optimisers.setup(Optimisers.AdamW(lr, (0.9, 0.999), 1e-5), ps)
    best_val = Inf
    best_ps = deepcopy(ps)
    hist_rows = NamedTuple[]

    for ep in 1:epochs
        # Cosine annealing schedule
        lr_ep = lr * (0.05 + 0.95 * 0.5 * (1.0 + cos(π * (ep - 1) / max(epochs, 1))))
        Optimisers.adjust!(opt, lr_ep)

        loss, grads = Zygote.withgradient(ps) do p
            t_hat, C_gl = pinn_c_eval_batch(p, tr)
            rel_err = (t_hat .- tr.t_true) ./ tr.t_true
            l_data = mean(abs2, rel_err) + 0.25 * mean(abs.(rel_err))

            # Physics monotonicity penalties: ∂C/∂U >= 0 and ∂C/∂HR <= 0
            out_u  = vec(mlp_forward(p, tr.X_u))
            C_u    = tr.C_base_u .* exp.(0.15 .* tanh.(out_u))
            l_mu   = mean(max.(0.0, C_gl .- C_u) ./ (C_gl .+ 1e-12))

            out_hr = vec(mlp_forward(p, tr.X_hr))
            C_hr   = tr.C_base_hr .* exp.(0.15 .* tanh.(out_hr))
            l_mhr  = mean(max.(0.0, C_hr .- C_gl) ./ (C_gl .+ 1e-12))

            l_reg = 5e-5 * _l2_params(p)
            return l_data + 0.5 * (l_mu + l_mhr) + l_reg
        end

        opt, ps = Optimisers.update(opt, ps, grads[1])

        if va !== nothing
            t_v, _ = pinn_c_eval_batch(ps, va)
            v_mape = mape_pct(t_v, va.t_true)
            if v_mape < best_val
                best_val = v_mape
                best_ps = deepcopy(ps)
            end
            if log_history && (ep % 50 == 0 || ep == 1 || ep == epochs)
                t_tr, _ = pinn_c_eval_batch(ps, tr)
                push!(hist_rows, (model="pinn_c", epoch=ep, train_loss=loss,
                                  train_metric=mape_pct(t_tr, tr.t_true), val_metric=v_mape))
            end
        elseif log_history && (ep % 50 == 0 || ep == 1 || ep == epochs)
            t_tr, _ = pinn_c_eval_batch(ps, tr)
            push!(hist_rows, (model="pinn_c", epoch=ep, train_loss=loss,
                              train_metric=mape_pct(t_tr, tr.t_true), val_metric=NaN))
        end
    end

    final_ps = (va !== nothing) ? best_ps : ps
    return PINNCSurrogate(final_ps; alpha_max=0.15), hist_rows
end

# =============================================================================
# 3. TRAIN MODEL 2: PINN-2Eq (Coupled 2-Eq + Neural Biot Ts + Energy Residual)
# =============================================================================

function prepare_pinn_2eq_tensors(d::Dict{String, Any}; d_cut_m::Float64=1.5e-5)
    N = length(d["t_evap_s"])
    M = N * N_TAU

    X        = zeros(Float64, 10, M)
    T_bar    = zeros(Float64, M)
    Ts_true  = zeros(Float64, M)
    dT_biot0 = zeros(Float64, M)
    T_inf    = zeros(Float64, M)
    Y_inf_c  = zeros(Float64, M)
    d_m      = zeros(Float64, M)
    rg       = zeros(Float64, M)
    Dab      = zeros(Float64, M)
    Sh0      = zeros(Float64, M)
    Nu0      = zeros(Float64, M)
    kg       = zeros(Float64, M)
    rl       = zeros(Float64, M)
    cpl      = zeros(Float64, M)
    m_eff    = zeros(Float64, M)
    md_true  = zeros(Float64, M)
    dTdt_true = zeros(Float64, M)
    is_t0    = zeros(Float64, M)
    s_norm   = zeros(Float64, N_TAU, N) # (d/d0)^2 along trajectory

    # Precompute Gauss-Legendre s quadrature interpolation weights (for t_evap consistency)
    # Using fixed indices across trajectory to integrate 1 / |dV/dt(tau)| * |dV/dtau| -> or simply
    # exact trajectory integral: t_evap = \int_0^1 (-dV/dtau / (-dV_dt_hat)) dtau!
    # Wait: even cleaner! By trapezoidal rule on tau_grid, V_true(tau_k) intervals have known ΔV_k!
    # t_evap_hat = sum_{k=1}^{N_TAU-1} (V_k - V_{k+1}) / (0.5 * (-dV_dt_hat[k] - dV_dt_hat[k+1]))!
    dV_seg = zeros(Float64, N_TAU - 1, N)

    for i in 1:N
        d0_mm = d["d0_m"][i] * 1e3
        U     = d["U_ms"][i]
        T0    = d["T0_C"][i]
        Tinf  = d["T_inf_C"][i]
        HR    = d["HR"][i]
        sc    = PINNSurrogates.compute_static_thermo_scales(Tinf, HR, T0)

        for k in 1:N_TAU
            col = (i - 1) * N_TAU + k
            Vk  = d["V_mm3"][k, i]
            Tbk = d["T_bar"][k, i]
            ph  = PINNSurrogates.pinn_2eq_precompute_physics(Vk, Tbk, U, Tinf, sc.Y_inf)

            X[:, col] = collect(pinn_2eq_features(
                ph.d_mm, d0_mm, Tbk, Tinf, HR, sc.T_wb, ph.Re, ph.B_M_Tb, ph.dT_biot0
            ))
            T_bar[col]     = Tbk
            Ts_true[col]   = d["Ts"][k, i]
            dT_biot0[col]  = ph.dT_biot0
            T_inf[col]     = Tinf
            Y_inf_c[col]   = ph.Y_inf_c
            d_m[col]       = ph.d_m
            rg[col]        = ph.rg
            Dab[col]       = ph.Dab
            Sh0[col]       = ph.Sh0
            Nu0[col]       = ph.Nu0
            kg[col]        = ph.kg
            rl[col]        = ph.rl
            cpl[col]       = ph.cpl
            d_eff          = sqrt(ph.d_m^2 + d_cut_m^2)
            m_eff[col]     = ph.rl * (π / 6.0) * (d_eff^3)
            md_true[col]   = d["m_dot"][k, i]
            dTdt_true[col] = d["dT_dt"][k, i]
            is_t0[col]     = (k == 1) ? 1.0 : 0.0
            s_norm[k, i]   = (ph.d_mm / d0_mm)^2
        end

        # Segment volumes [mm^3] between collocation points (including final drop to 0)
        for k in 1:(N_TAU - 1)
            dV_seg[k, i] = max(d["V_mm3"][k, i] - d["V_mm3"][k + 1, i], 0.0)
        end
    end

    return (
        N=N, n_tau=N_TAU, M=M, X=X, T_bar=T_bar, Ts_true=Ts_true,
        dT_biot0=dT_biot0, T_inf=T_inf, Y_inf_c=Y_inf_c,
        d_m=d_m, rg=rg, Dab=Dab, Sh0=Sh0, Nu0=Nu0, kg=kg, rl=rl, cpl=cpl,
        m_eff=m_eff, md_true=md_true, dTdt_true=dTdt_true, is_t0=is_t0,
        dV_seg=dV_seg, t_true=Float64.(d["t_evap_s"])
    )
end

function pinn_2eq_forward_batch(ps, tens)
    out = mlp_forward(ps, tens.X)
    o_Ts = view(out, 1, :)
    o_M  = view(out, 2, :)
    o_T  = view(out, 3, :)

    biot_fac = 0.05 .+ 1.20 .* PINNSurrogates.sigmoid_act.(o_Ts)
    Ts_hat = tens.T_bar .+ tens.dT_biot0 .* biot_fac
    Ts_K = Ts_hat .+ T_ZERO_C

    Y_s = surface_mass_fraction_from_psat.(psat_water.(Ts_K), P_ATM)
    B_M = max.((Y_s .- tens.Y_inf_c) ./ max.(1.0 .- Y_s, 1e-12), -0.999)
    DHev = latent_heat_water.(Ts_K)

    phi_M = exp.(0.40 .* tanh.(o_M))
    m_dot_hat = phi_M .* π .* tens.d_m .* tens.rg .* tens.Dab .* tens.Sh0 .* log1p.(B_M)
    dV_dt_hat = -(m_dot_hat ./ tens.rl) .* 1e9

    q_g = π .* tens.d_m .* tens.kg .* tens.Nu0 .* log1p_over_x.(B_M) .* (tens.T_inf .- Ts_hat)
    phi_T = exp.(0.50 .* tanh.(o_T))
    dT_dt_hat = phi_T .* (q_g .- m_dot_hat .* DHev) ./ (tens.m_eff .* tens.cpl)

    # Trajectory lifetime integral along V(tau): Δt_k ≈ ΔV_k / avg_rate_k
    rate_mat = reshape(-dV_dt_hat, tens.n_tau, tens.N)
    avg_rate = 0.5 .* (view(rate_mat, 1:(tens.n_tau - 1), :) .+ view(rate_mat, 2:tens.n_tau, :)) .+ 1e-12
    # Scale by 1 / (1 - tau_end_residual) ≈ 1/0.995 since TAU_GRID ends at 0.995
    t_hat = vec(sum(tens.dV_seg ./ avg_rate; dims=1)) ./ 0.995

    return (Ts_hat=Ts_hat, m_dot_hat=m_dot_hat, dT_dt_hat=dT_dt_hat, t_hat=t_hat)
end

function train_pinn_2eq(train_d, val_d; epochs::Int=700, lr::Float64=2.5e-3, seed::Int=42, log_history::Bool=false)
    rng = MersenneTwister(seed)
    ps = init_mlp_params(rng, [10, 32, 32, 3]; out_scale=0.01)
    tr = prepare_pinn_2eq_tensors(train_d)
    va = val_d === nothing ? nothing : prepare_pinn_2eq_tensors(val_d)

    opt = Optimisers.setup(Optimisers.AdamW(lr, (0.9, 0.999), 1e-5), ps)
    best_val = Inf
    best_ps = deepcopy(ps)
    hist_rows = NamedTuple[]

    w_t0 = 1.0 .+ 6.0 .* tr.is_t0

    for ep in 1:epochs
        lr_ep = lr * (0.05 + 0.95 * 0.5 * (1.0 + cos(π * (ep - 1) / max(epochs, 1))))
        Optimisers.adjust!(opt, lr_ep)

        loss, grads = Zygote.withgradient(ps) do p
            pr = pinn_2eq_forward_batch(p, tr)
            # 1. Surface temperature Biot closure loss (extra weight at t=0)
            l_Ts = mean(w_t0 .* abs2.(pr.Ts_hat .- tr.Ts_true))
            # 2. Mass flux relative loss (extra weight at t=0 for m_dot_0)
            rel_m = (pr.m_dot_hat .- tr.md_true) ./ (abs.(tr.md_true) .+ 1e-11)
            l_md = mean(w_t0 .* abs2.(rel_m))
            # 3. Energy balance derivative loss dT_bar/dt
            rel_dT = (pr.dT_dt_hat .- tr.dTdt_true) ./ (abs.(tr.dTdt_true) .+ 0.25)
            l_dT = mean(abs2, rel_dT)
            # 4. Lifetime integration consistency loss
            rel_t = (pr.t_hat .- tr.t_true) ./ tr.t_true
            l_tev = mean(abs2, rel_t) + 0.25 * mean(abs.(rel_t))

            l_reg = 1e-6 * _l2_params(p)
            return 0.5 * l_Ts + 1.5 * l_md + 0.4 * l_dT + 4.0 * l_tev + l_reg
        end

        opt, ps = Optimisers.update(opt, ps, grads[1])

        if va !== nothing
            pv = pinn_2eq_forward_batch(ps, va)
            v_score = mape_pct(pv.t_hat, va.t_true) + 0.25 * mape_pct(pv.m_dot_hat, va.md_true)
            if v_score < best_val
                best_val = v_score
                best_ps = deepcopy(ps)
            end
            if log_history && (ep % 50 == 0 || ep == 1 || ep == epochs)
                pt = pinn_2eq_forward_batch(ps, tr)
                push!(hist_rows, (model="pinn_2eq", epoch=ep, train_loss=loss,
                                  train_metric=mape_pct(pt.t_hat, tr.t_true),
                                  val_metric=mape_pct(pv.t_hat, va.t_true)))
            end
        elseif log_history && (ep % 50 == 0 || ep == 1 || ep == epochs)
            pt = pinn_2eq_forward_batch(ps, tr)
            push!(hist_rows, (model="pinn_2eq", epoch=ep, train_loss=loss,
                              train_metric=mape_pct(pt.t_hat, tr.t_true), val_metric=NaN))
        end
    end

    final_ps = (va !== nothing) ? best_ps : ps
    return PINN2EqSurrogate(final_ps; d_cut_m=1.5e-5), hist_rows
end

# =============================================================================
# 4. TRAIN MODEL 3: PI-DEEPONET (Solver-Free Physics-Informed Operator)
# =============================================================================

function prepare_deeponet_tensors(d::Dict{String, Any})
    N = length(d["t_evap_s"])
    taus = collect(TAU_GRID)
    X_cond = zeros(Float64, 10, N)
    t_base = zeros(Float64, N)
    T0_row = zeros(Float64, N)
    Twb_row = zeros(Float64, N)
    d0_mm_vec = zeros(Float64, N)

    # Physics inputs for initial flux m_dot_0 and Ts biot scaling along trajectory
    m0_phys_base = zeros(Float64, N)
    dT_biot0_mat = zeros(Float64, N, N_TAU)
    s_true       = zeros(Float64, N, N_TAU)
    Tb_true      = zeros(Float64, N, N_TAU)
    Ts_true      = zeros(Float64, N, N_TAU)
    dV_dt_true   = zeros(Float64, N, N_TAU)
    dT_dt_true   = zeros(Float64, N, N_TAU)

    for i in 1:N
        d0_mm = d["d0_m"][i] * 1e3
        U     = d["U_ms"][i]
        T0    = d["T0_C"][i]
        Tinf  = d["T_inf_C"][i]
        HR    = d["HR"][i]
        sc    = PINNSurrogates.compute_static_thermo_scales(Tinf, HR, T0)
        dT    = max(sc.dT, 1e-3)

        invC = sum(PINNSurrogates.GL_W[q] / calc_C_model7_smooth(d0_mm * PINNSurrogates.GL_SQRT_S[q], U, Tinf, HR, sc.T_wb, T0)
                   for q in eachindex(PINNSurrogates.GL_W))
        t_base[i]    = (sqrt(π) * (d0_mm^2) / (4.0 * dT)) * invC
        X_cond[:, i] = collect(deeponet_cond_features(d0_mm, U, T0, Tinf, HR, sc.T_wb, sc.Theta0, sc.BM_wb))
        T0_row[i]    = T0
        Twb_row[i]   = sc.T_wb
        d0_mm_vec[i] = d0_mm

        for k in 1:N_TAU
            Vk  = d["V_mm3"][k, i]
            Tbk = d["T_bar"][k, i]
            ph  = PINNSurrogates.pinn_2eq_precompute_physics(Vk, Tbk, U, Tinf, sc.Y_inf)
            dT_biot0_mat[i, k] = ph.dT_biot0
            if k == 1
                m0_phys_base[i] = π * ph.d_m * ph.rg * ph.Dab * ph.Sh0
            end
            s_true[i, k]     = (d["d_mm"][k, i] / d0_mm)^2
            Tb_true[i, k]    = Tbk
            Ts_true[i, k]    = d["Ts"][k, i]
            dV_dt_true[i, k] = d["dV_dt"][k, i]
            dT_dt_true[i, k] = d["dT_dt"][k, i]
        end
    end

    eps_tau = 1e-3
    T_tau   = hcat([collect(deeponet_trunk_features(t)) for t in taus]...)
    T_tau_p = hcat([collect(deeponet_trunk_features(min(t + eps_tau, 1.0))) for t in taus]...)
    T_tau_m = hcat([collect(deeponet_trunk_features(max(t - eps_tau, 0.0))) for t in taus]...)

    return (
        N=N, n_tau=N_TAU, P=16, eps_tau=eps_tau,
        taus=reshape(taus, 1, N_TAU),
        taus_p=reshape(clamp.(taus .+ eps_tau, 0.0, 1.0), 1, N_TAU),
        taus_m=reshape(clamp.(taus .- eps_tau, 0.0, 1.0), 1, N_TAU),
        T_tau=T_tau, T_tau_p=T_tau_p, T_tau_m=T_tau_m,
        X_cond=X_cond, t_base=t_base, t_true=Float64.(d["t_evap_s"]),
        T0=reshape(T0_row, N, 1), Twb=reshape(Twb_row, N, 1),
        V0=reshape((π / 6.0) .* (d0_mm_vec .^ 3), N, 1),
        dT_biot0=dT_biot0_mat, m0_phys_base=m0_phys_base,
        Y_inf=[PINNSurrogates.compute_static_thermo_scales(d["T_inf_C"][i], d["HR"][i], d["T0_C"][i]).Y_inf for i in 1:N],
        md0_true=Float64.(d["m_dot"][1, :]),
        s_true=s_true, Tb_true=Tb_true, Ts_true=Ts_true,
        dV_dt_true=dV_dt_true, dT_dt_true=dT_dt_true
    )
end

function deeponet_forward_batch(ps_b, ps_t, tens)
    P = tens.P
    B_all = mlp_forward(ps_b, tens.X_cond; act=PINNSurrogates.swish_act) # (3 + 3P, N)
    b_t   = vec(view(B_all, 1, :))
    b_lam = view(B_all, 2, :)
    b_md  = vec(view(B_all, 3, :))

    B_d  = view(B_all, 4:(3 + P), :)
    B_T  = view(B_all, (4 + P):(3 + 2P), :)
    B_Ts = view(B_all, (4 + 2P):(3 + 3P), :)

    Phi   = mlp_forward(ps_t, tens.T_tau;   act=tanh) # (P, n_tau)
    Phi_p = mlp_forward(ps_t, tens.T_tau_p; act=tanh)
    Phi_m = mlp_forward(ps_t, tens.T_tau_m; act=tanh)

    t_hat = tens.t_base .* exp.(0.18 .* tanh.(b_t))
    lam   = reshape(2.0 .+ exp.(2.5 .+ 1.5 .* tanh.(b_lam)), tens.N, 1)

    tau   = tens.taus
    tau_p = tens.taus_p
    tau_m = tens.taus_m

    # Hard BCs trajectories
    dot_d  = B_d' * Phi
    dot_dp = B_d' * Phi_p
    dot_dm = B_d' * Phi_m
    s_hat  = (1.0 .- tau) .+ (tau .* (1.0 .- tau)) .* dot_d
    s_p    = (1.0 .- tau_p) .+ (tau_p .* (1.0 .- tau_p)) .* dot_dp
    s_m    = (1.0 .- tau_m) .+ (tau_m .* (1.0 .- tau_m)) .* dot_dm

    dot_T  = B_T' * Phi
    dot_Tp = B_T' * Phi_p
    dot_Tm = B_T' * Phi_m
    Tb_hat = tens.Twb .+ (tens.T0 .- tens.Twb) .* exp.(-lam .* tau) .+ tau .* dot_T
    Tb_p   = tens.Twb .+ (tens.T0 .- tens.Twb) .* exp.(-lam .* tau_p) .+ tau_p .* dot_Tp
    Tb_m   = tens.Twb .+ (tens.T0 .- tens.Twb) .* exp.(-lam .* tau_m) .+ tau_m .* dot_Tm

    dot_Ts = B_Ts' * Phi
    Ts_hat = Tb_hat .+ tens.dT_biot0 .* (0.05 .+ 1.20 .* PINNSurrogates.sigmoid_act.(dot_Ts))

    # Physics collocation time derivatives: dV/dt = V0 * (3/2) * sqrt(s) * (ds/dtau) / t_hat
    dtau_denom = max.(tau_p .- tau_m, 1e-6)
    ds_dtau = (s_p .- s_m) ./ dtau_denom
    dT_dtau = (Tb_p .- Tb_m) ./ dtau_denom
    t_mat   = reshape(t_hat, tens.N, 1)
    dV_dt_pred = tens.V0 .* (1.5 .* sqrt.(max.(s_hat, 1e-4)) .* ds_dtau) ./ t_mat
    dT_dt_pred = dT_dtau ./ t_mat

    # Initial flux m_dot_0 at tau=0
    Ts0_hat = view(Ts_hat, :, 1)
    Y_s0    = surface_mass_fraction_from_psat.(psat_water.(Ts0_hat .+ T_ZERO_C), P_ATM)
    B_M0    = max.((Y_s0 .- tens.Y_inf) ./ max.(1.0 .- Y_s0, 1e-12), -0.999)
    md0_hat = exp.(0.40 .* tanh.(b_md)) .* tens.m0_phys_base .* log1p.(B_M0)

    return (
        t_hat=t_hat, s_hat=s_hat, Tb_hat=Tb_hat, Ts_hat=Ts_hat,
        dV_dt_pred=dV_dt_pred, dT_dt_pred=dT_dt_pred, md0_hat=md0_hat
    )
end

function train_pi_deeponet(train_d, val_d; epochs::Int=800, lr::Float64=2.5e-3, seed::Int=42, log_history::Bool=false)
    rng = MersenneTwister(seed)
    P = 16
    ps_b = init_mlp_params(rng, [10, 48, 48, 3 + 3 * P]; out_scale=0.02)
    ps_t = init_mlp_params(rng, [4, 24, 24, P]; out_scale=0.05)
    ps = (branch=ps_b, trunk=ps_t)

    tr = prepare_deeponet_tensors(train_d)
    va = val_d === nothing ? nothing : prepare_deeponet_tensors(val_d)

    opt = Optimisers.setup(Optimisers.AdamW(lr, (0.9, 0.999), 1e-5), ps)
    best_val = Inf
    best_ps = deepcopy(ps)
    hist_rows = NamedTuple[]

    for ep in 1:epochs
        lr_ep = lr * (0.05 + 0.95 * 0.5 * (1.0 + cos(π * (ep - 1) / max(epochs, 1))))
        Optimisers.adjust!(opt, lr_ep)

        loss, grads = Zygote.withgradient(ps) do p
            pr = deeponet_forward_batch(p.branch, p.trunk, tr)
            # 1. Lifetime loss
            rel_t = (pr.t_hat .- tr.t_true) ./ tr.t_true
            l_tev = mean(abs2, rel_t) + 0.25 * mean(abs.(rel_t))
            # 2. Trajectory state matching (d^2/d0^2, T_bar, Ts)
            l_s  = mean(abs2, pr.s_hat .- tr.s_true)
            l_Tb = mean(abs2, (pr.Tb_hat .- tr.Tb_true) ./ 3.0)
            l_Ts = mean(abs2, (pr.Ts_hat .- tr.Ts_true) ./ 3.0)
            # 3. Initial mass flux m_dot_0
            l_md0 = mean(abs2, (pr.md0_hat .- tr.md0_true) ./ (abs.(tr.md0_true) .+ 1e-11))
            # 4. Physics ODE collocation residual on interior tau (indices 2..n_tau-1)
            idx_in = 2:(tr.n_tau - 1)
            v_scale = abs.(view(tr.dV_dt_true, :, idx_in)) .+ 1e-6
            t_scale = abs.(view(tr.dT_dt_true, :, idx_in)) .+ 0.25
            l_ode_v = mean(abs2, (view(pr.dV_dt_pred, :, idx_in) .- view(tr.dV_dt_true, :, idx_in)) ./ v_scale)
            l_ode_T = mean(abs2, (view(pr.dT_dt_pred, :, idx_in) .- view(tr.dT_dt_true, :, idx_in)) ./ t_scale)

            l_reg = 5e-5 * (_l2_params(p.branch) + _l2_params(p.trunk))
            return 4.0 * l_tev + 2.0 * l_s + 2.5 * l_Tb + 1.5 * l_Ts + 1.5 * l_md0 + 0.4 * l_ode_v + 0.3 * l_ode_T + l_reg
        end

        opt, ps = Optimisers.update(opt, ps, grads[1])

        if va !== nothing
            pv = deeponet_forward_batch(ps.branch, ps.trunk, va)
            v_mape = mape_pct(pv.t_hat, va.t_true)
            if v_mape < best_val
                best_val = v_mape
                best_ps = deepcopy(ps)
            end
            if log_history && (ep % 50 == 0 || ep == 1 || ep == epochs)
                pt = deeponet_forward_batch(ps.branch, ps.trunk, tr)
                push!(hist_rows, (model="pi_deeponet", epoch=ep, train_loss=loss,
                                  train_metric=mape_pct(pt.t_hat, tr.t_true), val_metric=v_mape))
            end
        elseif log_history && (ep % 50 == 0 || ep == 1 || ep == epochs)
            pt = deeponet_forward_batch(ps.branch, ps.trunk, tr)
            push!(hist_rows, (model="pi_deeponet", epoch=ep, train_loss=loss,
                              train_metric=mape_pct(pt.t_hat, tr.t_true), val_metric=NaN))
        end
    end

    final_ps = (va !== nothing) ? best_ps : ps
    return PIDeepONetSurrogate(final_ps.branch, final_ps.trunk; p_basis=P, alpha_t=0.18), hist_rows
end

# =============================================================================
# 5. TRAIN MODEL 4: RESMLP (Pure Data-Driven Baseline — Ablation)
# =============================================================================

function prepare_resmlp_tensors(d::Dict{String, Any}; x_mean=nothing, x_std=nothing)
    N = length(d["t_evap_s"])
    X_raw = zeros(Float64, 5, N)
    for i in 1:N
        X_raw[1, i] = log(max(d["d0_m"][i] * 1e3, 1e-4))
        X_raw[2, i] = d["U_ms"][i]
        X_raw[3, i] = d["T0_C"][i]
        X_raw[4, i] = d["T_inf_C"][i]
        X_raw[5, i] = d["HR"][i]
    end
    xm = x_mean === nothing ? vec(mean(X_raw; dims=2)) : x_mean
    xs = x_std  === nothing ? max.(vec(std(X_raw; dims=2)), 1e-6) : x_std
    X = (X_raw .- xm) ./ xs
    return (
        X=X, x_mean=xm, x_std=xs,
        y_logt = log.(max.(Float64.(d["t_evap_s"]), 1e-6)),
        t_true = Float64.(d["t_evap_s"]),
        y_init = vcat(reshape((Float64.(d["Ts"][1, :]) .- 25.0) ./ 15.0, 1, N),
                      reshape(log.(max.(Float64.(d["m_dot"][1, :]), 1e-20)), 1, N))
    )
end

function train_resmlp(train_d, val_d; epochs::Int=900, lr::Float64=3e-3, seed::Int=42, log_history::Bool=false)
    rng = MersenneTwister(seed)
    tr = prepare_resmlp_tensors(train_d)
    va = val_d === nothing ? nothing : prepare_resmlp_tensors(val_d; x_mean=tr.x_mean, x_std=tr.x_std)

    ps_tevap = init_mlp_params(rng, [5, 64, 64, 64, 1]; out_scale=0.1)
    ps_init  = init_mlp_params(rng, [5, 48, 48, 2]; out_scale=0.1)
    ps = (tevap=ps_tevap, init=ps_init)

    opt = Optimisers.setup(Optimisers.AdamW(lr, (0.9, 0.999), 1e-5), ps)
    best_val = Inf
    best_ps = deepcopy(ps)
    hist_rows = NamedTuple[]

    for ep in 1:epochs
        lr_ep = lr * (0.05 + 0.95 * 0.5 * (1.0 + cos(π * (ep - 1) / max(epochs, 1))))
        Optimisers.adjust!(opt, lr_ep)

        loss, grads = Zygote.withgradient(ps) do p
            logt_hat = vec(mlp_forward(p.tevap, tr.X))
            l_t = mean(abs2, logt_hat .- tr.y_logt)
            init_hat = mlp_forward(p.init, tr.X)
            l_i = mean(abs2, init_hat .- tr.y_init)
            l_reg = 1e-5 * (_l2_params(p.tevap) + _l2_params(p.init))
            return l_t + 0.5 * l_i + l_reg
        end

        opt, ps = Optimisers.update(opt, ps, grads[1])

        if va !== nothing
            t_v = exp.(clamp.(vec(mlp_forward(ps.tevap, va.X)), -10.0, 12.0))
            v_mape = mape_pct(t_v, va.t_true)
            if v_mape < best_val
                best_val = v_mape
                best_ps = deepcopy(ps)
            end
            if log_history && (ep % 50 == 0 || ep == 1 || ep == epochs)
                t_tr = exp.(clamp.(vec(mlp_forward(ps.tevap, tr.X)), -10.0, 12.0))
                push!(hist_rows, (model="resmlp", epoch=ep, train_loss=loss,
                                  train_metric=mape_pct(t_tr, tr.t_true), val_metric=v_mape))
            end
        elseif log_history && (ep % 50 == 0 || ep == 1 || ep == epochs)
            t_tr = exp.(clamp.(vec(mlp_forward(ps.tevap, tr.X)), -10.0, 12.0))
            push!(hist_rows, (model="resmlp", epoch=ep, train_loss=loss,
                              train_metric=mape_pct(t_tr, tr.t_true), val_metric=NaN))
        end
    end

    final_ps = (va !== nothing) ? best_ps : ps
    return ResMLPSurrogate(final_ps.tevap, final_ps.init, tr.x_mean, tr.x_std), hist_rows
end

# =============================================================================
# 6. EVALUATION HELPERS (ACTUAL ODE SOLVE / OPERATOR INFERENCE)
# =============================================================================

function evaluate_dataset_metrics(model_id::Symbol, model_obj, d::Dict{String, Any})
    N = length(d["t_evap_s"])
    t_hat   = Vector{Float64}(undef, N)
    Ts0_hat = Vector{Float64}(undef, N)
    md0_hat = Vector{Float64}(undef, N)
    d_mae   = Float64[]
    Tb_mae  = Float64[]

    for i in 1:N
        d0_mm = d["d0_m"][i] * 1e3
        U     = d["U_ms"][i]
        T0    = d["T0_C"][i]
        Tinf  = d["T_inf_C"][i]
        HR    = d["HR"][i]
        sc    = PINNSurrogates.compute_static_thermo_scales(Tinf, HR, T0)

        if model_id === :pinn_c
            sol = solve_pinn_c(model_obj, d0_mm, U, T0, Tinf, HR)
            t_hat[i] = Float64(sol.t[end])
            Ts0_hat[i] = sc.T_wb
            C0 = pinn_c_eval_C(model_obj, d0_mm, d0_mm, U, T0, Tinf, HR, sc.T_wb, sc.Theta0, sc.BM_wb)
            dV0 = -C0 * sqrt(sphere_surface_from_d(d0_mm)) * (Tinf - sc.T_wb)
            md0_hat[i] = -dV0 * 1e-9 * liquid_density_water(sc.T_wb + T_ZERO_C)
            # Sample diameter MAE on trajectory
            d_preds = [sphere_diameter_from_volume(max(sol(clamp(TAU_GRID[k] * t_hat[i], 0.0, t_hat[i]))[1], 0.0)) for k in 1:N_TAU]
            push!(d_mae, mean(abs.(d_preds .- d["d_mm"][:, i])))
        elseif model_id === :pinn_2eq
            sol = solve_pinn_2eq(model_obj, d0_mm, U, T0, Tinf, HR)
            t_hat[i] = Float64(sol.t[end])
            V0 = (π / 6.0) * (d0_mm^3)
            st0 = pinn_2eq_step_physics(model_obj, V0, T0, d0_mm, U, Tinf, HR, sc.T_wb, sc.Y_inf)
            Ts0_hat[i] = Float64(st0.Ts_hat)
            md0_hat[i] = Float64(st0.m_dot)
            # Trajectory MAE for both d_mm and T_bar
            d_preds = Float64[]
            T_preds = Float64[]
            for k in 1:N_TAU
                tk = clamp(TAU_GRID[k] * d["t_evap_s"][i], 0.0, t_hat[i])
                uk = sol(tk)
                push!(d_preds, sphere_diameter_from_volume(max(uk[1], 0.0)))
                push!(T_preds, uk[2])
            end
            push!(d_mae, mean(abs.(d_preds .- d["d_mm"][:, i])))
            push!(Tb_mae, mean(abs.(T_preds .- d["T_bar"][:, i])))
        elseif model_id === :pi_deeponet
            op = deeponet_eval_case(model_obj, d0_mm, U, T0, Tinf, HR; taus=TAU_GRID)
            t_hat[i] = op.t_evap_s
            Ts0_hat[i] = op.T_s0_C
            md0_hat[i] = op.m_dot_0_kgs
            push!(d_mae, mean(abs.(op.d_mm .- d["d_mm"][:, i])))
            push!(Tb_mae, mean(abs.(op.T_bar .- d["T_bar"][:, i])))
        elseif model_id === :resmlp
            op = resmlp_eval_case(model_obj, d0_mm, U, T0, Tinf, HR)
            t_hat[i] = op.t_evap_s
            Ts0_hat[i] = op.T_s0_C
            md0_hat[i] = op.m_dot_0_kgs
        end
    end

    t_true   = Float64.(d["t_evap_s"])
    Ts0_true = Float64.(d["Ts"][1, :])
    md0_true = Float64.(d["m_dot"][1, :])

    return (
        t_evap_MAPE_pct = mape_pct(t_hat, t_true),
        t_evap_MaxARE_pct = max_are_pct(t_hat, t_true),
        Ts0_MAE_C = mae_val(Ts0_hat, Ts0_true),
        mdot0_MAPE_pct = mape_pct(md0_hat, md0_true),
        traj_d_MAE_mm = isempty(d_mae) ? NaN : mean(d_mae),
        traj_Tbar_MAE_C = isempty(Tb_mae) ? NaN : mean(Tb_mae),
    )
end

# =============================================================================
# 7. MAIN WORKFLOW: 5-FOLD CV ON V2 + FULL FIT + V1 OOD EVALUATION
# =============================================================================

function main()
    args = parse_cli()
    smoke = args["smoke"]
    regen = args["regen-data"]

    println("="^80)
    println("JAIIO 2026 PINN & DEEPONET TRAINING SUITE (Julia + Zygote + AdamW)")
    println("Timestamp: ", Dates.now(), " | Smoke: ", smoke)
    println("="^80)

    v2_conds = get_v2_conditions()
    v1_conds = get_v1_conditions()
    if smoke
        v2_conds = v2_conds[1:21:end][1:min(120, end)]
        v1_conds = v1_conds[1:10:end][1:min(60, end)]
    end

    v2_cache = joinpath(DATA_DIR, smoke ? "pinn_v2_traj_smoke.npz" : "pinn_v2_trajectories.npz")
    v1_cache = joinpath(DATA_DIR, smoke ? "pinn_v1_traj_smoke.npz" : "pinn_v1_trajectories.npz")

    v2_d = build_or_load_trajectories("V2", v2_conds, v2_cache; force_regen=regen)
    v1_d = build_or_load_trajectories("V1", v1_conds, v1_cache; force_regen=regen)

    N_v2 = length(v2_d["t_evap_s"])
    rng = MersenneTwister(2026)
    perm = randperm(rng, N_v2)

    ep_c    = smoke ? 120 : 850
    ep_2eq  = smoke ? 120 : 650
    ep_deep = smoke ? 120 : 750
    ep_res  = smoke ? 150 : 900

    # --- A. 5-Fold Cross-Validation on V2 ---
    k_folds = 5
    cv_rows = NamedTuple[]
    println("\n--- Part A: 5-Fold Cross-Validation across V2 (N=$N_v2 cases) ---")

    for fold in 1:k_folds
        val_idx = perm[fold:k_folds:N_v2]
        trn_idx = setdiff(perm, val_idx)
        trn_d = slice_dataset(v2_d, trn_idx)
        val_d = slice_dataset(v2_d, val_idx)
        @printf("\nFold %d/%d (train=%d, val=%d):\n", fold, k_folds, length(trn_idx), length(val_idx))

        m_c, _    = train_pinn_c(trn_d, val_d; epochs=ep_c, seed=100 + fold)
        met_c     = evaluate_dataset_metrics(:pinn_c, m_c, val_d)
        @printf("  PINN-1Eq    | t_evap MAPE=%.3f%% | mdot0 MAPE=%.1f%% | d MAE=%.4f mm\n",
                met_c.t_evap_MAPE_pct, met_c.mdot0_MAPE_pct, met_c.traj_d_MAE_mm)
        push!(cv_rows, merge((model="pinn_c", fold=fold, n=length(val_idx)), met_c))

        m_2eq, _  = train_pinn_2eq(trn_d, val_d; epochs=ep_2eq, seed=200 + fold)
        met_2eq   = evaluate_dataset_metrics(:pinn_2eq, m_2eq, val_d)
        @printf("  PINN-2Eq    | t_evap MAPE=%.3f%% | mdot0 MAPE=%.2f%% | Ts0 MAE=%.3f°C | Tbar MAE=%.3f°C\n",
                met_2eq.t_evap_MAPE_pct, met_2eq.mdot0_MAPE_pct, met_2eq.Ts0_MAE_C, met_2eq.traj_Tbar_MAE_C)
        push!(cv_rows, merge((model="pinn_2eq", fold=fold, n=length(val_idx)), met_2eq))

        m_deep, _ = train_pi_deeponet(trn_d, val_d; epochs=ep_deep, seed=300 + fold)
        met_deep  = evaluate_dataset_metrics(:pi_deeponet, m_deep, val_d)
        @printf("  PI-DeepONet | t_evap MAPE=%.3f%% | mdot0 MAPE=%.2f%% | Ts0 MAE=%.3f°C | Tbar MAE=%.3f°C\n",
                met_deep.t_evap_MAPE_pct, met_deep.mdot0_MAPE_pct, met_deep.Ts0_MAE_C, met_deep.traj_Tbar_MAE_C)
        push!(cv_rows, merge((model="pi_deeponet", fold=fold, n=length(val_idx)), met_deep))

        m_res, _  = train_resmlp(trn_d, val_d; epochs=ep_res, seed=400 + fold)
        met_res   = evaluate_dataset_metrics(:resmlp, m_res, val_d)
        @printf("  ResMLP-Puro | t_evap MAPE=%.3f%% | mdot0 MAPE=%.2f%% | Ts0 MAE=%.3f°C\n",
                met_res.t_evap_MAPE_pct, met_res.mdot0_MAPE_pct, met_res.Ts0_MAE_C)
        push!(cv_rows, merge((model="resmlp", fold=fold, n=length(val_idx)), met_res))
    end

    df_cv = DataFrame(cv_rows)
    CSV.write(joinpath(RESULTS_DIR, "pinn_cv_metrics.csv"), df_cv)

    df_cv_sum = combine(groupby(df_cv, :model),
        :t_evap_MAPE_pct   => mean => :cv_tevap_MAPE_pct,
        :t_evap_MaxARE_pct => maximum => :cv_tevap_MaxARE_pct,
        :Ts0_MAE_C         => mean => :cv_Ts0_MAE_C,
        :mdot0_MAPE_pct    => mean => :cv_mdot0_MAPE_pct,
        :traj_d_MAE_mm     => mean => :cv_traj_d_MAE_mm,
        :traj_Tbar_MAE_C   => mean => :cv_traj_Tbar_MAE_C,
    )
    println("\n--- 5-Fold CV Mean Summary (V2) ---")
    println(df_cv_sum)

    # --- B. Train Final Models on Full V2 Dataset (with 10% internal monitor for history curves) ---
    println("\n--- Part B: Final Production Training on V2 (with V1 Out-of-Distribution Test) ---")
    # For clean convergence plots, log on 90/10 split first, then train on 100% V2 (or use full V2 with V1 monitor!)
    # We use 95/5 split solely to avoid overfitting; V1 remains 100% unseen OOD!
    n_val_mon = max(20, N_v2 ÷ 10)
    val_mon_idx = perm[1:n_val_mon]
    trn_mon_idx = perm[(n_val_mon + 1):N_v2]
    trn_mon = slice_dataset(v2_d, trn_mon_idx)
    val_mon = slice_dataset(v2_d, val_mon_idx)

    final_c, h_c       = train_pinn_c(trn_mon, val_mon; epochs=ep_c, seed=2026, log_history=true)
    final_2eq, h_2eq   = train_pinn_2eq(trn_mon, val_mon; epochs=ep_2eq, seed=2027, log_history=true)
    final_deep, h_deep = train_pi_deeponet(trn_mon, val_mon; epochs=ep_deep, seed=2028, log_history=true)
    final_res, h_res   = train_resmlp(trn_mon, val_mon; epochs=ep_res, seed=2029, log_history=true)

    bundle = PINNSurrogateBundle(final_c, final_2eq, final_deep, final_res)
    save_pinn_bundle(bundle, PINN_WEIGHTS_PATH)
    println("Saved trained PINN weights to: ", PINN_WEIGHTS_PATH)

    df_hist = DataFrame(vcat(h_c, h_2eq, h_deep, h_res))
    CSV.write(joinpath(RESULTS_DIR, "pinn_training_history.csv"), df_hist)

    # --- C. Out-of-Distribution Generalization on V1 (640 unseen cases, up to Tinf=65°C, HR=10-90%) ---
    ood_rows = NamedTuple[]
    for (m_id, m_obj) in ((:pinn_c, final_c), (:pinn_2eq, final_2eq), (:pi_deeponet, final_deep), (:resmlp, final_res))
        mv2 = evaluate_dataset_metrics(m_id, m_obj, v2_d)
        mv1 = evaluate_dataset_metrics(m_id, m_obj, v1_d)
        push!(ood_rows, (
            model = string(m_id),
            display_name = display_name(m_id),
            v2_tevap_MAPE_pct = mv2.t_evap_MAPE_pct,
            v1_ood_tevap_MAPE_pct = mv1.t_evap_MAPE_pct,
            v1_ood_tevap_MaxARE_pct = mv1.t_evap_MaxARE_pct,
            v1_ood_Ts0_MAE_C = mv1.Ts0_MAE_C,
            v1_ood_mdot0_MAPE_pct = mv1.mdot0_MAPE_pct,
            v1_ood_traj_d_MAE_mm = mv1.traj_d_MAE_mm,
            v1_ood_traj_Tbar_MAE_C = mv1.traj_Tbar_MAE_C,
        ))
    end
    df_ood = DataFrame(ood_rows)
    CSV.write(joinpath(RESULTS_DIR, "pinn_v1_ood_metrics.csv"), df_ood)
    println("\n--- Out-of-Distribution Generalization (V1: 640 cases) ---")
    println(df_ood)
end

main()
