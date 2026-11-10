"""
    Evaporacion.Neural.PINNSurrogates

Physics-Informed Neural Network (PINN) and Operator Learning surrogates for
spherical droplet evaporation in Julia:
  1. PINNCSurrogate      (:pinn_c)      - 1-Eq Physics-Informed Residual d²-Law ODE
  2. PINN2EqSurrogate    (:pinn_2eq)    - 2-Eq Coupled [V(t), T_bar(t)] + Neural Biot Ts Closure
  3. PIDeepONetSurrogate (:pi_deeponet) - Solver-Free Physics-Informed DeepONet Operator
  4. ResMLPSurrogate     (:resmlp)      - Pure Data-Driven Residual MLP (Ablation Baseline)

Zero data leakage: all models depend strictly on physical inputs (d0, U, T0, T_inf, HR)
and never receive t_evap as an input feature.
"""
module PINNSurrogates

using LinearAlgebra
using NPZ
using OrdinaryDiffEq
using QuadGK
using Random
using StaticArrays

using ...Common
using ...Droplet
using ...Surrogates

export StaticMLP, to_static_mlp, mlp_forward, static_mlp_eval, init_mlp_params,
       PINNCSurrogate, PINN2EqSurrogate, PIDeepONetSurrogate, ResMLPSurrogate,
       PINNSurrogateBundle,
       calc_C_model7_smooth,
       pinn_c_features, pinn_c_eval_C, pinn_c_predict_tevap, solve_pinn_c,
       pinn_2eq_features, pinn_2eq_step_physics, solve_pinn_2eq,
       deeponet_cond_features, deeponet_trunk_features, deeponet_eval_case,
       resmlp_eval_case,
       save_pinn_bundle, load_pinn_bundle, get_default_pinn_bundle

# =============================================================================
# ACTIVATIONS & STATIC/VECTORIZED MLP KERNELS
# =============================================================================

@inline swish_act(x::T) where {T<:Real} = x / (one(T) + exp(-x))
@inline sigmoid_act(x::T) where {T<:Real} = one(T) / (one(T) + exp(-x))

const _GL_N_RAW, _GL_W_RAW = QuadGK.gauss(8)
const GL_S = 0.5 .* (_GL_N_RAW .+ 1.0)
const GL_W = 0.5 .* _GL_W_RAW
const GL_SQRT_S = sqrt.(GL_S)

"""
    init_mlp_params(rng, dims::Vector{Int}; out_scale=0.05)

Initializes a pure-functional NamedTuple parameter structure `(layers = ((W=..., b=...), ...))`
suitable for Zygote + Optimisers.jl training in Float64.
"""
function init_mlp_params(rng::AbstractRNG, dims::Vector{Int}; out_scale::Float64=0.05)
    n_layers = length(dims) - 1
    layers = ntuple(n_layers) do i
        d_in, d_out = dims[i], dims[i + 1]
        scale = (i == n_layers) ? out_scale : sqrt(2.0 / (d_in + d_out))
        W = randn(rng, Float64, d_out, d_in) .* scale
        b = zeros(Float64, d_out)
        (W=W, b=b)
    end
    return (layers=layers,)
end

"""
    mlp_forward(ps, X::AbstractMatrix; act=swish_act)

Vectorized forward pass on feature matrix X of shape (d_in, N_batch).
"""
function mlp_forward(ps::NamedTuple, X::AbstractMatrix; act=swish_act)
    h = X
    n = length(ps.layers)
    for i in 1:(n - 1)
        W = ps.layers[i].W
        b = ps.layers[i].b
        h = act.(W * h .+ b)
    end
    W_last = ps.layers[n].W
    b_last = ps.layers[n].b
    return W_last * h .+ b_last
end

struct StaticMLP{L<:Tuple, F}
    layers::L
    act::F
end

function to_static_mlp(ps::NamedTuple; act=swish_act)
    slayers = map(ps.layers) do l
        d_out, d_in = size(l.W)
        (W = SMatrix{d_out, d_in, Float64}(l.W), b = SVector{d_out, Float64}(l.b))
    end
    return StaticMLP(slayers, act)
end

@inline function _eval_static_layers(layers::Tuple{T1}, h, act) where {T1}
    W, b = layers[1].W, layers[1].b
    return W * h + b
end

@inline function _eval_static_layers(layers::Tuple{T1, T2, Vararg}, h, act) where {T1, T2}
    W, b = layers[1].W, layers[1].b
    h1 = act.(W * h + b)
    return _eval_static_layers(Base.tail(layers), h1, act)
end

@inline function static_mlp_eval(m::StaticMLP, x::SVector)
    return _eval_static_layers(m.layers, x, m.act)
end

# =============================================================================
# SMOOTH C∞ PHYSICS BASELINE (NO PIECEWISE REYNOLDS DISCONTINUITIES)
# =============================================================================

"""
    calc_C_model7_smooth(d_mm, U, T_inf, HR, T_wb, T0)

Smooth C∞-continuous blend of the calibrated `modelo_7` regimes across Re=10 and Re=100,
eliminating piecewise jumps as the physics prior for PINN-C and PI-DeepONet.
"""
@inline function calc_C_model7_smooth(d_mm::T, U::Real, T_inf::Real, HR::Real, T_wb::Real, T0::Real) where {T<:Real}
    Re = calc_reynolds_surrogate(d_mm, U, T_inf)
    # Smooth rational/sigmoid regime gates around Re=10 and Re=100
    w_mid  = sigmoid_act((Re - T(10.0)) / T(2.5))
    w_high = sigmoid_act((Re - T(100.0)) / T(15.0))
    w1 = one(T) - w_mid
    w2 = w_mid * (one(T) - w_high)
    w3 = w_high

    t_ratio = max(T(T_inf) / T(25.0), T(1e-6))
    dT = T(T_inf) - T(T_wb)
    Theta = abs(dT) > T(1e-3) ? (T(T0) - T(T_wb)) / dT : zero(T)
    hr_clamped = clamp(T(HR), zero(T), T(0.999))
    one_m_hr = one(T) - hr_clamped
    du = d_mm * abs(T(U))

    # Calibrated V2 modelo_7 regime values
    c_low  = begin
        A0, n0, al = T(5.346202e-05), T(0.281287), T(-0.050473)
        A1, n1, be = T(9.376442e-05), T(0.167952), T(-0.028636)
        ga, de, ci = T(0.008057), T(0.000816), T(0.495)
        (A0 * (t_ratio^n0) * (one_m_hr^al) + A1 * (t_ratio^n1) * (du^ci) * (one_m_hr^be)) * exp(ga * Theta + de * Theta^2)
    end
    c_mid  = begin
        A0, n0, al = T(5.292208e-05), T(0.281515), T(-0.054463)
        A1, n1, be = T(9.533896e-05), T(0.166045), T(-0.027342)
        ga, de, ci = T(0.007001), T(0.000643), T(0.496)
        (A0 * (t_ratio^n0) * (one_m_hr^al) + A1 * (t_ratio^n1) * (du^ci) * (one_m_hr^be)) * exp(ga * Theta + de * Theta^2)
    end
    c_high = begin
        A0, n0, al = T(5.145208e-05), T(0.283085), T(-0.056330)
        A1, n1, be = T(9.684243e-05), T(0.167299), T(-0.023637)
        ga, de, ci = T(0.006643), T(0.000690), T(0.493)
        (A0 * (t_ratio^n0) * (one_m_hr^al) + A1 * (t_ratio^n1) * (du^ci) * (one_m_hr^be)) * exp(ga * Theta + de * Theta^2)
    end

    return w1 * c_low + w2 * c_mid + w3 * c_high
end

# Helper to compute thermodynamic scale inputs shared across models
@inline function compute_static_thermo_scales(T_inf::Real, HR::Real, T0::Real)
    T_wb = stull_twb(Float64(T_inf), Float64(HR) * 100.0)
    dT = Float64(T_inf) - T_wb
    Theta0 = abs(dT) > 1e-3 ? clamp((Float64(T0) - T_wb) / dT, -5.0, 5.0) : 0.0
    psat_wb = psat_water(T_wb + T_ZERO_C)
    Ys_wb = surface_mass_fraction_from_psat(psat_wb, P_ATM)
    psat_inf = psat_water(Float64(T_inf) + T_ZERO_C)
    Y_inf = surface_mass_fraction_from_psat(clamp(Float64(HR), 0.0, 1.0) * psat_inf, P_ATM)
    BM_wb = max((Ys_wb - Y_inf) / max(1.0 - Ys_wb, 1e-12), 0.0)
    return (T_wb=T_wb, dT=dT, Theta0=Theta0, Y_inf=Y_inf, BM_wb=BM_wb)
end

# =============================================================================
# 1. PINN-C: 1-EQUATION PHYSICS-INFORMED RESIDUAL d²-LAW ODE
# =============================================================================

struct PINNCSurrogate{M<:StaticMLP}
    ps::NamedTuple
    smlp::M
    alpha_max::Float64
end

function PINNCSurrogate(ps::NamedTuple; alpha_max::Float64=0.45)
    return PINNCSurrogate(ps, to_static_mlp(ps; act=swish_act), alpha_max)
end

@inline function pinn_c_features(
    d_mm::T, d0_mm::Real, U::Real, T0::Real, T_inf::Real, HR::Real,
    Theta0::Real, BM_wb::Real
) where {T<:Real}
    Re = calc_reynolds_surrogate(d_mm, U, T_inf)
    f1 = log1p(max(Re, zero(T))) / T(5.0)
    f2 = Re / (T(10.0) + Re)
    f3 = Re / (T(100.0) + Re)
    f4 = (T(T_inf) - T(30.0)) / T(20.0)
    f5 = (log(max(one(T) - T(HR), T(1e-4))) + T(0.7)) / T(0.5)
    f6 = T(Theta0) / T(2.0)
    f7 = log1p(max(T(BM_wb), zero(T))) * T(25.0)
    f8 = d_mm / max(T(d0_mm), T(1e-4)) - T(0.5)
    f9 = (log(max(T(d0_mm), T(0.05))) + T(0.5)) / T(1.5)
    return SVector{9, T}(f1, f2, f3, f4, f5, f6, f7, f8, f9)
end

@inline function pinn_c_eval_C(
    model::PINNCSurrogate, d_mm::T, d0_mm::Real, U::Real, T0::Real, T_inf::Real,
    HR::Real, T_wb::Real, Theta0::Real, BM_wb::Real
) where {T<:Real}
    C_base = calc_C_model7_smooth(d_mm, U, T_inf, HR, T_wb, T0)
    x = pinn_c_features(d_mm, d0_mm, U, T0, T_inf, HR, Theta0, BM_wb)
    delta = static_mlp_eval(model.smlp, x)[1]
    factor = exp(T(model.alpha_max) * tanh(delta))
    return C_base * factor
end

"""
    pinn_c_predict_tevap(model::PINNCSurrogate, d0_mm, U, T0, T_inf, HR)

Analytically integrates t_evap via 8-point Gauss-Legendre quadrature (zero ODE stepping).
"""
function pinn_c_predict_tevap(model::PINNCSurrogate, d0_mm::Real, U::Real, T0::Real, T_inf::Real, HR::Real)
    sc = compute_static_thermo_scales(T_inf, HR, T0)
    dT = max(sc.dT, 1e-3)
    invC = 0.0
    @inbounds for q in eachindex(GL_W)
        dq = d0_mm * GL_SQRT_S[q]
        Cq = pinn_c_eval_C(model, dq, d0_mm, U, T0, T_inf, HR, sc.T_wb, sc.Theta0, sc.BM_wb)
        invC += GL_W[q] / Cq
    end
    return (sqrt(π) * (Float64(d0_mm)^2) / (4.0 * dT)) * invC
end

struct PINNCRHSParams{M<:PINNCSurrogate}
    model::M
    d0_mm::Float64
    U::Float64
    T0::Float64
    T_inf::Float64
    HR::Float64
    T_wb::Float64
    Theta0::Float64
    BM_wb::Float64
end

function pinn_c_ode_rhs!(du, u, p::PINNCRHSParams, t)
    V = u[1]
    if V <= 1e-9
        du[1] = 0.0
        return nothing
    end
    d_mm = sphere_diameter_from_volume(V)
    C = pinn_c_eval_C(p.model, d_mm, p.d0_mm, p.U, p.T0, p.T_inf, p.HR, p.T_wb, p.Theta0, p.BM_wb)
    S_mm2 = sphere_surface_from_d(d_mm)
    dT = p.T_inf - p.T_wb
    du[1] = -C * sqrt(S_mm2) * dT
    return nothing
end

function solve_pinn_c(
    model::PINNCSurrogate, d0_mm::Float64, U::Float64, T0::Float64, T_inf::Float64, HR::Float64;
    t_end::Float64=30000.0, solver=Tsit5(), rtol::Float64=1e-6, atol::Float64=1e-9
)
    sc = compute_static_thermo_scales(T_inf, HR, T0)
    p = PINNCRHSParams(model, d0_mm, U, T0, T_inf, HR, sc.T_wb, sc.Theta0, sc.BM_wb)
    V0 = (π / 6.0) * (d0_mm^3)
    cb = ContinuousCallback((u, t, integrator) -> u[1] - 1e-9, terminate!)
    prob = ODEProblem(pinn_c_ode_rhs!, [V0], (0.0, t_end), p)
    return solve(prob, solver; callback=cb, reltol=rtol, abstol=atol)
end

# =============================================================================
# 2. PINN-2Eq: COUPLED 2-EQUATION THERMO-DIFFUSIVE PINN ([V(t), T_bar(t)] + Ts)
# =============================================================================

struct PINN2EqSurrogate{M<:StaticMLP}
    ps::NamedTuple
    smlp::M
    d_cut_m::Float64
end

function PINN2EqSurrogate(ps::NamedTuple; d_cut_m::Float64=1.5e-5)
    return PINN2EqSurrogate(ps, to_static_mlp(ps; act=swish_act), d_cut_m)
end

"""
    pinn_2eq_precompute_physics(V_mm3, T_bar, U, T_inf, Y_inf)

Computes instantaneous gas/liquid film properties and the 1st-order analytical Biot-Clausius
surface temperature shift ΔT_biot0 (zero root-finding iterations).
"""
@inline function pinn_2eq_precompute_physics(
    V_mm3::T, T_bar::Real, U::Real, T_inf::Real, Y_inf::Real
) where {T<:Real}
    d_mm = sphere_diameter_from_volume(max(V_mm3, T(1e-12)))
    d_m = d_mm * T(1e-3)
    T_bar_C = T(T_bar)
    T_bar_K = T_bar_C + T(T_ZERO_C)
    T_inf_C = T(T_inf)
    T_film_K = T(0.5) * (T_inf_C + T_bar_C) + T(T_ZERO_C)

    # Liquid properties at T_bar
    rl  = liquid_density_water(T_bar_K)
    cpl = liquid_cp_water(T_bar_K)
    kl  = liquid_thermal_conductivity_water(T_bar_K)

    # Gas film properties at 1/2 film rule
    rg  = gas_density_air(T_film_K, P_ATM)
    mug = gas_viscosity_air(T_film_K)
    cpg = gas_cp_air(T_film_K)
    kg  = gas_thermal_conductivity_air(T_film_K)
    Dab = gas_diffusivity_water_air(T_film_K, P_ATM)

    Re  = calc_re(d_m, U, rg, mug)
    Pr  = calc_pr(mug, cpg, kg)
    Sc  = calc_sc(mug, rg, Dab)
    Nu0 = nusselt_froessling(Re, Pr)
    Sh0 = sherwood_froessling(Re, Sc)

    P_sat_Tb = psat_water(T_bar_K)
    Y_s_Tb   = surface_mass_fraction_from_psat(P_sat_Tb, P_ATM)
    Y_inf_c  = clamp(T(Y_inf), zero(T), T(0.99999))
    B_M_Tb   = max((Y_s_Tb - Y_inf_c) / max(one(T) - Y_s_Tb, T(1e-12)), T(-0.999))
    DHev_Tb  = latent_heat_water(T_bar_K)

    # Scale-invariant energy residual at Ts = T_bar divided by 5*kl:
    # R_0 = 0.5*Nu0*kg*(T_inf - T_bar) - 0.5*Sh0*rg*Dab*ln(1+BM_Tb)*DHev_Tb
    R0 = T(0.5) * Nu0 * kg * (T_inf_C - T_bar_C) - T(0.5) * Sh0 * rg * Dab * log1p(B_M_Tb) * DHev_Tb
    dT_biot0 = R0 / (T(5.0) * kl)

    return (
        d_mm=d_mm, d_m=d_m, rl=rl, cpl=cpl, kl=kl,
        rg=rg, kg=kg, Dab=Dab, Re=Re, Nu0=Nu0, Sh0=Sh0,
        Y_inf_c=Y_inf_c, B_M_Tb=B_M_Tb, dT_biot0=dT_biot0
    )
end

@inline function pinn_2eq_features(
    d_mm::T, d0_mm::Real, T_bar::Real, T_inf::Real, HR::Real, T_wb::Real,
    Re::T, B_M_Tb::T, dT_biot0::T
) where {T<:Real}
    f1 = log1p(max(Re, zero(T))) / T(5.0)
    f2 = Re / (T(10.0) + Re)
    f3 = Re / (T(100.0) + Re)
    f4 = (T(T_inf) - T(30.0)) / T(20.0)
    f5 = (T(T_bar) - T(25.0)) / T(20.0)
    f6 = (T(T_bar) - T(T_wb)) / T(20.0)
    f7 = (log(max(one(T) - T(HR), T(1e-4))) + T(0.7)) / T(0.5)
    f8 = log1p(max(B_M_Tb, T(-0.5))) * T(20.0)
    f9 = d_mm / max(T(d0_mm), T(1e-4)) - T(0.5)
    f10 = clamp(dT_biot0 / T(5.0), T(-2.0), T(2.0))
    return SVector{10, T}(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10)
end

"""
    pinn_2eq_step_physics(model::PINN2EqSurrogate, V_mm3, T_bar, d0_mm, U, T_inf, HR, T_wb, Y_inf)

Evaluates the coupled physics + neural closure at state (V_mm3, T_bar):
Returns `(dV_dt, dT_dt, Ts_hat, m_dot)`.
"""
@inline function pinn_2eq_step_physics(
    model::PINN2EqSurrogate,
    V_mm3::T, T_bar::Real, d0_mm::Real, U::Real, T_inf::Real, HR::Real, T_wb::Real, Y_inf::Real
) where {T<:Real}
    ph = pinn_2eq_precompute_physics(V_mm3, T_bar, U, T_inf, Y_inf)
    x = pinn_2eq_features(ph.d_mm, d0_mm, T_bar, T_inf, HR, T_wb, ph.Re, ph.B_M_Tb, ph.dT_biot0)
    out = static_mlp_eval(model.smlp, x)

    o_Ts = out[1]
    o_M  = out[2]
    o_T  = out[3]

    # 1. Surface temperature closure: Ts_hat = T_bar + dT_biot0 * scale(o_Ts)
    # Since internal circulation (bq >= 5) & Clausius-Clapeyron slope shrink |Ts - T_bar|, scale in (0.05, 1.25)
    biot_factor = T(0.05) + T(1.20) * sigmoid_act(o_Ts)
    Ts_hat = T(T_bar) + ph.dT_biot0 * biot_factor

    # 2. Mass transfer at Ts_hat with neural blowing/film multiplier
    Ts_K = Ts_hat + T(T_ZERO_C)
    P_sat_s = psat_water(Ts_K)
    Y_s = surface_mass_fraction_from_psat(P_sat_s, P_ATM)
    B_M = max((Y_s - ph.Y_inf_c) / max(one(T) - Y_s, T(1e-12)), T(-0.999))
    DHev = latent_heat_water(Ts_K)

    phi_M = exp(T(0.40) * tanh(o_M))
    m_dot = phi_M * T(π) * ph.d_m * ph.rg * ph.Dab * ph.Sh0 * log1p(B_M)
    dV_dt = -(m_dot / ph.rl) * T(1e9) # mm^3/s

    # 3. Energy balance equation with softly regularized thermal mass denominator
    q_g = T(π) * ph.d_m * ph.kg * ph.Nu0 * log1p_over_x(B_M) * (T(T_inf) - Ts_hat)
    d_eff = sqrt(ph.d_m^2 + T(model.d_cut_m)^2)
    m_eff = ph.rl * (T(π) / T(6.0)) * (d_eff^3)
    phi_T = exp(T(0.50) * tanh(o_T))
    dT_dt = phi_T * (q_g - m_dot * DHev) / (m_eff * ph.cpl)

    return (dV_dt=dV_dt, dT_dt=dT_dt, Ts_hat=Ts_hat, m_dot=m_dot)
end

struct PINN2EqRHSParams{M<:PINN2EqSurrogate}
    model::M
    d0_mm::Float64
    U::Float64
    T_inf::Float64
    HR::Float64
    T_wb::Float64
    Y_inf::Float64
end

function pinn_2eq_ode_rhs!(du, u, p::PINN2EqRHSParams, t)
    V = u[1]
    T_bar = u[2]
    if V <= 1e-9
        du[1] = 0.0
        du[2] = 0.0
        return nothing
    end
    st = pinn_2eq_step_physics(p.model, V, T_bar, p.d0_mm, p.U, p.T_inf, p.HR, p.T_wb, p.Y_inf)
    du[1] = st.dV_dt
    du[2] = st.dT_dt
    return nothing
end

function solve_pinn_2eq(
    model::PINN2EqSurrogate, d0_mm::Float64, U::Float64, T0::Float64, T_inf::Float64, HR::Float64;
    t_end::Float64=30000.0, solver=AutoTsit5(Rosenbrock23(autodiff=false)), rtol::Float64=1e-6, atol::Float64=1e-9
)
    sc = compute_static_thermo_scales(T_inf, HR, T0)
    p = PINN2EqRHSParams(model, d0_mm, U, T_inf, HR, sc.T_wb, sc.Y_inf)
    V0 = (π / 6.0) * (d0_mm^3)
    cb = ContinuousCallback((u, t, integrator) -> u[1] - 1e-9, terminate!)
    prob = ODEProblem(pinn_2eq_ode_rhs!, [V0, T0], (0.0, t_end), p)
    return solve(prob, solver; callback=cb, reltol=rtol, abstol=[atol, atol])
end

# =============================================================================
# 3. PI-DEEPONET: SOLVER-FREE PHYSICS-INFORMED OPERATOR NETWORK
# =============================================================================

struct PIDeepONetSurrogate{B<:StaticMLP, T<:StaticMLP}
    ps_branch::NamedTuple
    ps_trunk::NamedTuple
    sbranch::B
    strunk::T
    p_basis::Int
    alpha_t::Float64
end

function PIDeepONetSurrogate(ps_branch::NamedTuple, ps_trunk::NamedTuple; p_basis::Int=16, alpha_t::Float64=0.40)
    return PIDeepONetSurrogate(
        ps_branch, ps_trunk,
        to_static_mlp(ps_branch; act=swish_act),
        to_static_mlp(ps_trunk; act=tanh),
        p_basis, alpha_t
    )
end

@inline function deeponet_cond_features(
    d0_mm::Real, U::Real, T0::Real, T_inf::Real, HR::Real,
    T_wb::Real, Theta0::Real, BM_wb::Real
)
    Re0 = calc_reynolds_surrogate(Float64(d0_mm), U, T_inf)
    return SVector{10, Float64}(
        (log(max(Float64(d0_mm), 0.05)) + 0.5) / 1.5,
        log1p(Re0) / 5.0,
        Re0 / (10.0 + Re0),
        Re0 / (100.0 + Re0),
        (Float64(T_inf) - 30.0) / 20.0,
        (Float64(T0) - 25.0) / 20.0,
        (log(max(1.0 - Float64(HR), 1e-4)) + 0.7) / 0.5,
        Float64(Theta0) / 2.0,
        log1p(max(Float64(BM_wb), 0.0)) * 25.0,
        (Float64(T_inf) - Float64(T_wb)) / 20.0
    )
end

@inline function deeponet_trunk_features(tau::T) where {T<:Real}
    tc = clamp(tau, zero(T), one(T))
    return SVector{4, T}(tc, sqrt(tc + T(1e-8)), log1p(T(25.0) * tc) / T(3.26), one(T) - tc)
end

"""
    deeponet_eval_case(model::PIDeepONetSurrogate, d0_mm, U, T0, T_inf, HR; taus=[0.0, ...])

Single-pass solver-free operator evaluation (`nfev = 1`):
Returns `(t_evap_s, taus, d_mm, T_bar, Ts, m_dot_0, dV_dt_bulk)`.
"""
function deeponet_eval_case(
    model::PIDeepONetSurrogate, d0_mm::Float64, U::Float64, T0::Float64, T_inf::Float64, HR::Float64;
    taus::AbstractVector{Float64}=Float64[0.0, 0.01, 0.03, 0.08, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95, 1.0]
)
    sc = compute_static_thermo_scales(T_inf, HR, T0)
    dT = max(sc.dT, 1e-3)

    # Base physics t_evap from smooth modelo_7 via 8-pt GL (still nfev=1 analytic evaluation)
    invC = 0.0
    @inbounds for q in eachindex(GL_W)
        invC += GL_W[q] / calc_C_model7_smooth(d0_mm * GL_SQRT_S[q], U, T_inf, HR, sc.T_wb, T0)
    end
    t_base = (sqrt(π) * (d0_mm^2) / (4.0 * dT)) * invC

    cond_x = deeponet_cond_features(d0_mm, U, T0, T_inf, HR, sc.T_wb, sc.Theta0, sc.BM_wb)
    b = static_mlp_eval(model.sbranch, cond_x)

    P = model.p_basis
    b_t   = b[1]
    b_lam = b[2]
    b_md  = b[3]

    t_evap = t_base * exp(model.alpha_t * tanh(b_t))
    lam = 2.0 + exp(2.5 + 1.5 * tanh(b_lam))

    n_tau = length(taus)
    d_arr  = Vector{Float64}(undef, n_tau)
    Tb_arr = Vector{Float64}(undef, n_tau)
    Ts_arr = Vector{Float64}(undef, n_tau)

    @inbounds for i in 1:n_tau
        tau = clamp(taus[i], 0.0, 1.0)
        phi = static_mlp_eval(model.strunk, deeponet_trunk_features(tau))

        # Inner products over basis slices
        dot_d  = 0.0
        dot_T  = 0.0
        dot_Ts = 0.0
        for j in 1:P
            phij = phi[j]
            dot_d  += b[3 + j] * phij
            dot_T  += b[3 + P + j] * phij
            dot_Ts += b[3 + 2 * P + j] * phij
        end

        # Hard BCs: d(0)=d0, d(1)=0; T(0)=T0
        s_tau = clamp((1.0 - tau) + tau * (1.0 - tau) * dot_d, 0.0, 1.4)
        d_i = d0_mm * sqrt(s_tau)
        Tb_i = sc.T_wb + (T0 - sc.T_wb) * exp(-lam * tau) + tau * dot_T

        # Biot Ts closure anchored to instantaneous physical R0 at (d_i, Tb_i)
        V_i = (π / 6.0) * (max(d_i, 1e-4)^3)
        ph_i = pinn_2eq_precompute_physics(V_i, Tb_i, U, T_inf, sc.Y_inf)
        Ts_i = Tb_i + ph_i.dT_biot0 * (0.05 + 1.20 * sigmoid_act(dot_Ts))

        d_arr[i]  = d_i
        Tb_arr[i] = Tb_i
        Ts_arr[i] = Ts_i
    end

    # Initial flux m_dot_0 from physical Sherwood law evaluated at Ts_hat(0)
    V0 = (π / 6.0) * (d0_mm^3)
    ph0 = pinn_2eq_precompute_physics(V0, T0, U, T_inf, sc.Y_inf)
    # Evaluate Ts at tau=0
    phi0 = static_mlp_eval(model.strunk, deeponet_trunk_features(0.0))
    dot_Ts0 = sum(b[3 + 2 * P + j] * phi0[j] for j in 1:P)
    Ts0 = T0 + ph0.dT_biot0 * (0.05 + 1.20 * sigmoid_act(dot_Ts0))

    Ts0_K = Ts0 + T_ZERO_C
    Y_s0 = surface_mass_fraction_from_psat(psat_water(Ts0_K), P_ATM)
    B_M0 = max((Y_s0 - ph0.Y_inf_c) / max(1.0 - Y_s0, 1e-12), -0.999)
    m_dot_0 = exp(0.40 * tanh(b_md)) * π * ph0.d_m * ph0.rg * ph0.Dab * ph0.Sh0 * log1p(B_M0)
    dV_dt_bulk = V0 / max(t_evap, 1e-12)

    return (
        t_evap_s = t_evap,
        taus = collect(taus),
        d_mm = d_arr,
        T_bar = Tb_arr,
        Ts = Ts_arr,
        T_s0_C = Ts0,
        m_dot_0_kgs = m_dot_0,
        dV_dt_bulk_mm3s = dV_dt_bulk,
        T_wb_C = sc.T_wb
    )
end

# =============================================================================
# 4. RESMLP: PURE DATA-DRIVEN BASELINE (ABLATION — NO PHYSICS PRIOR/RESIDUALS)
# =============================================================================

struct ResMLPSurrogate{M1<:StaticMLP, M2<:StaticMLP}
    ps_tevap::NamedTuple
    ps_init::NamedTuple
    s_tevap::M1
    s_init::M2
    x_mean::SVector{5, Float64}
    x_std::SVector{5, Float64}
end

function ResMLPSurrogate(ps_tevap::NamedTuple, ps_init::NamedTuple, x_mean::AbstractVector, x_std::AbstractVector)
    return ResMLPSurrogate(
        ps_tevap, ps_init,
        to_static_mlp(ps_tevap; act=swish_act),
        to_static_mlp(ps_init; act=swish_act),
        SVector{5, Float64}(x_mean), SVector{5, Float64}(x_std)
    )
end

function resmlp_eval_case(model::ResMLPSurrogate, d0_mm::Float64, U::Float64, T0::Float64, T_inf::Float64, HR::Float64)
    raw = SVector{5, Float64}(log(max(d0_mm, 1e-4)), U, T0, T_inf, HR)
    xn = (raw .- model.x_mean) ./ model.x_std
    log_t = static_mlp_eval(model.s_tevap, xn)[1]
    t_evap = exp(clamp(log_t, -10.0, 12.0))
    init_out = static_mlp_eval(model.s_init, xn)
    Ts0 = init_out[1] * 15.0 + 25.0
    m_dot_0 = exp(clamp(init_out[2], -35.0, -5.0))
    V0 = (π / 6.0) * (d0_mm^3)
    return (
        t_evap_s = t_evap,
        T_s0_C = Ts0,
        m_dot_0_kgs = m_dot_0,
        dV_dt_bulk_mm3s = V0 / max(t_evap, 1e-12),
    )
end

# =============================================================================
# SERIALIZATION (NPZ) & DEFAULT CACHE
# =============================================================================

struct PINNSurrogateBundle
    pinn_c::PINNCSurrogate
    pinn_2eq::PINN2EqSurrogate
    pi_deeponet::PIDeepONetSurrogate
    resmlp::ResMLPSurrogate
end

function _pack_mlp!(dict::Dict{String, Any}, prefix::String, ps::NamedTuple)
    dict["$(prefix).n_layers"] = Int64[length(ps.layers)]
    for (i, l) in enumerate(ps.layers)
        dict["$(prefix).W$(i)"] = Array(l.W)
        dict["$(prefix).b$(i)"] = Array(l.b)
    end
end

function _unpack_mlp(dict::AbstractDict{String}, prefix::String)
    n_layers = Int(dict["$(prefix).n_layers"][1])
    layers = ntuple(n_layers) do i
        (W = Float64.(dict["$(prefix).W$(i)"]), b = Float64.(dict["$(prefix).b$(i)"]))
    end
    return (layers=layers,)
end

function save_pinn_bundle(bundle::PINNSurrogateBundle, path::AbstractString)
    dict = Dict{String, Any}()
    _pack_mlp!(dict, "pinn_c", bundle.pinn_c.ps)
    dict["pinn_c.alpha_max"] = Float64[bundle.pinn_c.alpha_max]

    _pack_mlp!(dict, "pinn_2eq", bundle.pinn_2eq.ps)
    dict["pinn_2eq.d_cut_m"] = Float64[bundle.pinn_2eq.d_cut_m]

    _pack_mlp!(dict, "deeponet.branch", bundle.pi_deeponet.ps_branch)
    _pack_mlp!(dict, "deeponet.trunk", bundle.pi_deeponet.ps_trunk)
    dict["deeponet.meta"] = Float64[bundle.pi_deeponet.p_basis, bundle.pi_deeponet.alpha_t]

    _pack_mlp!(dict, "resmlp.tevap", bundle.resmlp.ps_tevap)
    _pack_mlp!(dict, "resmlp.init", bundle.resmlp.ps_init)
    dict["resmlp.x_mean"] = collect(bundle.resmlp.x_mean)
    dict["resmlp.x_std"]  = collect(bundle.resmlp.x_std)

    mkpath(dirname(path))
    tmp_path = path * ".tmp.npz"
    npzwrite(tmp_path, dict)
    mv(tmp_path, path; force=true)
    return path
end

function load_pinn_bundle(path::AbstractString)
    dict = npzread(path)
    pinn_c = PINNCSurrogate(
        _unpack_mlp(dict, "pinn_c");
        alpha_max = Float64(dict["pinn_c.alpha_max"][1])
    )
    pinn_2eq = PINN2EqSurrogate(
        _unpack_mlp(dict, "pinn_2eq");
        d_cut_m = Float64(dict["pinn_2eq.d_cut_m"][1])
    )
    meta = dict["deeponet.meta"]
    pi_deeponet = PIDeepONetSurrogate(
        _unpack_mlp(dict, "deeponet.branch"),
        _unpack_mlp(dict, "deeponet.trunk");
        p_basis = Int(round(meta[1])), alpha_t = Float64(meta[2])
    )
    resmlp = ResMLPSurrogate(
        _unpack_mlp(dict, "resmlp.tevap"),
        _unpack_mlp(dict, "resmlp.init"),
        Float64.(dict["resmlp.x_mean"]),
        Float64.(dict["resmlp.x_std"])
    )
    return PINNSurrogateBundle(pinn_c, pinn_2eq, pi_deeponet, resmlp)
end

const _DEFAULT_BUNDLE_CACHE = Ref{Union{Nothing, PINNSurrogateBundle}}(nothing)

function get_default_pinn_bundle(path::AbstractString)
    if _DEFAULT_BUNDLE_CACHE[] !== nothing
        return _DEFAULT_BUNDLE_CACHE[]
    end
    isfile(path) || error("Trained PINN weights not found at $path. Run JAIIO/train_jaiio_pinns.jl first.")
    b = load_pinn_bundle(path)
    _DEFAULT_BUNDLE_CACHE[] = b
    return b
end

end # module PINNSurrogates
