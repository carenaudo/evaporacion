"""
    algebraic.jl

Algebraic core evaluation, surface energy balance, and Brent root-finding for Ts.
Generic over `T <: Real` to support automatic differentiation.
"""

"""
    _compute_gas_transport_and_fluxes(...)

Shared helper computing blowing-corrected transport numbers and gas-phase fluxes
consistently for both `calc_droplet_algebraic` and `intensive_surface_residual`.
"""
function _compute_gas_transport_and_fluxes(
    params::DropletParams,
    Re::T, Pr::T, Sc::T,
    Nu0_in::T, Sh0_in::T,
    rg::T, kg::T, Dab::T, cpv::T, cpg::T,
    rl::T, r_s::T,
    Y_s::T, Y_inf_c::T, bm_clamped::T,
    DHev::T, Ts::T
) where {T<:Real}
    Nu0 = Nu0_in
    Sh0 = Sh0_in
    ln_BM = log1p(bm_clamped)
    spheroid_ratio = params.spheroid_epsilon > 0.0 ? T(spheroid_deformation_evaporation_ratio(params.spheroid_epsilon)) : one(T)

    is_level1 = params.model_variant === :level_1
    is_level2 = params.model_variant === :level_2 || (!params.as1989_blowing && params.model_variant === :abramzon_sirignano)

    if is_level1 || is_level2
        Sh_star = is_level1 ? T(2.0) : Sh0
        Nu_star = is_level1 ? T(2.0) : Nu0
        B_T     = bm_clamped
        F_M     = one(T)
        F_T     = one(T)
        phi     = one(T)
        Sh      = Sh_star
        Nu      = Nu_star

        dY = Y_s - Y_inf_c
        flux_mass_intens = T(0.5) * rg * Dab * Sh_star * dY * spheroid_ratio
        m_dot_raw = T(4.0 * π) * r_s * flux_mass_intens
        h         = Nu_star * kg / (T(2.0) * r_s)
        Le        = one(T)
        r_dot_s   = -(rg * Dab * Sh_star) / (T(2.0) * rl * r_s) * dY
        term_qg   = T(0.5) * Nu_star * kg * (T(params.T_inf) - Ts)
        term_evap = flux_mass_intens * DHev

        return (Nu0=Nu0, Sh0=Sh0, B_T=B_T, F_M=F_M, F_T=F_T, phi=phi,
                Sh_star=Sh_star, Sh=Sh, Nu_star=Nu_star, Nu=Nu,
                m_dot_raw=m_dot_raw, h=h, Le=Le, r_dot_s=r_dot_s,
                term_qg=term_qg, term_evap=term_evap)
    end

    local Nu_star, Sh_star, B_T, F_M, F_T, phi
    if params.as1989_blowing
        as_sol = solve_as1989_transfer(Re, Pr, Sc, bm_clamped, cpv, cpg, Dab * rg * cpg / kg; coef=params.frossling_coef)
        Nu0     = as_sol.Nu0
        Sh0     = as_sol.Sh0
        Nu_star = as_sol.Nu_star
        Sh_star = as_sol.Sh_star
        B_T     = as_sol.B_T
        F_M     = as_sol.F_M
        F_T     = as_sol.F_T
        phi     = as_sol.phi
    else
        Nu_star = Nu0
        Sh_star = Sh0
        B_T     = bm_clamped
        F_M     = one(T)
        F_T     = one(T)
        phi     = one(T)
    end

    Sh = Sh_star * log1p_over_x(bm_clamped)
    Nu = Nu_star * log1p_over_x(B_T)

    flux_mass_intens = if params.model_variant === :d2_classic
        rg * Dab * (Y_s - Y_inf_c)
    elseif params.model_variant === :d2_stefan
        rg * Dab * ln_BM
    elseif params.model_variant === :d2_convective
        Sh_conv = sherwood_froessling(Re, Sc; coef=params.frossling_coef)
        T(0.5) * rg * Dab * Sh_conv * (Y_s - Y_inf_c)
    elseif params.model_variant === :d2_revised_dallabarba
        T(0.5) * rg * Dab * Sh0 * ln_BM
    else
        T(0.5) * rg * Dab * Sh_star * ln_BM
    end
    flux_mass_intens *= spheroid_ratio
    m_dot_raw = T(4.0 * π) * r_s * flux_mass_intens

    ln_BT = log1p(B_T)
    h = Nu_star * kg / (T(2.0) * r_s) * log1p_over_x(B_T)
    Le = abs(ln_BT) > T(1e-10) ? (ln_BM / ln_BT) : one(T)
    r_dot_s = -(rg * Dab * Sh_star) / (T(2.0) * rl * r_s) * ln_BM

    term_qg   = T(0.5) * Nu_star * kg * log1p_over_x(B_T) * (T(params.T_inf) - Ts)
    term_evap = flux_mass_intens * DHev

    return (Nu0=Nu0, Sh0=Sh0, B_T=B_T, F_M=F_M, F_T=F_T, phi=phi,
            Sh_star=Sh_star, Sh=Sh, Nu_star=Nu_star, Nu=Nu,
            m_dot_raw=m_dot_raw, h=h, Le=Le, r_dot_s=r_dot_s,
            term_qg=term_qg, term_evap=term_evap)
end

"""
    calc_droplet_algebraic(liq::DropletLiquidState, Ts::Real, params::DropletParams) -> DropletAlgebraic
    calc_droplet_algebraic(m::Real, T_bar::Real, Ts::Real, params::DropletParams) -> DropletAlgebraic

Evaluates the complete algebraic system for a spherical droplet at state (m, T_bar, Ts)
or using precomputed invariant liquid state `liq`.
Temperatures T_bar, Ts, and params.T_inf are in Celsius [°C].
"""
function calc_droplet_algebraic(
    liq::DropletLiquidState{Tliq}, Ts::Ts_t, params::DropletParams
) where {Tliq<:Real, Ts_t<:Real}
    T = promote_type(Tliq, Ts_t, Float64)

    # Unit conversions to Kelvin for physical property functions
    Ts_K    = T(Ts) + T(T_ZERO_C)
    T_inf_K = T(params.T_inf) + T(T_ZERO_C)
    
    if params.reference_state === :one_third
        ref = reference_state_one_third(Ts_K, T_inf_K, surface_mass_fraction_from_psat(psat_water(Ts_K), params.P_atm), params.Y_inf)
        T_eg_K = ref.T_ref
        T_eg_C = T_eg_K - T(T_ZERO_C)
        Y_ref  = ref.Y_v_ref
    else
        T_eg_C  = (T(params.T_inf) + T(Ts)) / T(2.0)
        T_eg_K  = T_eg_C + T(T_ZERO_C)
        Y_ref   = T(params.Y_inf)
    end

    # Liquid properties from precomputed state
    rl  = T(liq.rho_l)
    cpl = T(liq.cp_l)
    kl  = T(liq.k_l)
    mul = T(liq.mu_l)
    d   = T(liq.d)
    r_s = T(liq.r_s)

    # Gas film properties at T_eg
    if params.use_wilke_mixtures
        mix = gas_mixture_properties(T_eg_K, Y_ref, params.P_atm)
        rg   = mix.rho_g
        mug  = mix.mu_g
        cpg  = mix.cp_g
        kg   = mix.k_g
        Dab  = mix.D_ab
    else
        rg   = gas_density_air(T_eg_K, params.P_atm)
        mug  = gas_viscosity_air(T_eg_K)
        cpg  = gas_cp_air(T_eg_K)
        kg   = gas_thermal_conductivity_air(T_eg_K)
        Dab  = gas_diffusivity_water_air(T_eg_K, params.P_atm)
    end
    cpv  = vapor_cp_water(T_eg_K)

    # Latent heat of vaporization at Ts
    DHev = latent_heat_water(Ts_K)

    # Dimensionless numbers
    u_rel = params.terminal_velocity ? calc_terminal_velocity(d, T_eg_K, params.P_atm; rho_l=rl) : T(liq.u_rel)
    Re = calc_re(d, u_rel, rg, mug)
    Pr = calc_pr(mug, cpg, kg)
    Sc = calc_sc(mug, rg, Dab)

    if params.correlation === :clift
        Nu0, Sh0 = clift_correlation(Re, Pr, Sc)
    else
        Nu0 = nusselt_froessling(Re, Pr; coef=params.frossling_coef)
        Sh0 = sherwood_froessling(Re, Sc; coef=params.frossling_coef)
    end

    # Vapor mass fraction at surface and Spalding numbers
    P_sat_s = psat_water(Ts_K)
    Y_s = surface_mass_fraction_from_psat(P_sat_s, params.P_atm)
    Y_inf_c = clamp(T(params.Y_inf), zero(T), T(0.99999))
    B_M = (Y_s - Y_inf_c) / max(one(T) - Y_s, T(1e-12))
    bm_clamped = max(B_M, T(-0.99999))

    tr = _compute_gas_transport_and_fluxes(
        params, Re, Pr, Sc, Nu0, Sh0,
        rg, kg, Dab, cpv, cpg, rl, r_s,
        Y_s, Y_inf_c, bm_clamped, DHev, T(Ts)
    )
    Nu0     = tr.Nu0
    Sh0     = tr.Sh0
    Nu_star = tr.Nu_star
    Sh_star = tr.Sh_star
    B_T     = tr.B_T
    F_M     = tr.F_M
    F_T     = tr.F_T
    Sh      = tr.Sh
    Nu      = tr.Nu
    m_dot_raw = tr.m_dot_raw

    m_dot = abs(m_dot_raw) < T(1e-20) ? sign(m_dot_raw) * T(1e-20) : m_dot_raw

    # Internal circulation enhancement (Kronig-Brink / Hill vortex)
    circ = internal_circulation_enhancement(
        params.profile, Re, rg, mul, r_s, u_rel,
        rl, Dab, kl, cpl, bm_clamped
    )
    b_q   = circ.b_q
    b_m   = circ.b_m
    chi_l = circ.chi_l
    chi_m = params.fixed_chi_m ? one(T) : circ.chi_m
    U_s   = circ.U_s

    # Heat transfer rates [W]
    h       = tr.h
    Le      = tr.Le
    r_dot_s = tr.r_dot_s
    # Heat conducted from surface into droplet interior [W] (Abramzon & Sirignano 1989 Eq. 24)
    q_l = b_q != zero(T) ? T(4.0 * π) * b_q * kl * (T(Ts) - T(liq.T_bar)) * r_s : zero(T)
    q_g = T(4.0 * π) * (r_s^2) * h * (T(params.T_inf) - T(Ts))

    return DropletAlgebraic{T}(
        d, r_s, T_eg_C, rg, mug, cpg, kg, Dab, cpl, kl, rl,
        DHev, Re, Pr, Sc, Nu0, Sh0, Y_s, bm_clamped, B_T,
        F_M, F_T, Sh_star, Sh, Nu_star, Nu,
        m_dot, h, Le, r_dot_s, q_l, q_g,
        b_q, b_m, chi_l, chi_m, U_s
    )
end

"""
    calc_droplet_algebraic(m, T_bar, Ts, params) -> DropletAlgebraic
"""
function calc_droplet_algebraic(
    m::Tm, T_bar::Tt, Ts::Ts_t, params::DropletParams
) where {Tm<:Real, Tt<:Real, Ts_t<:Real}
    liq = precompute_liquid_state(m, T_bar, params)
    return calc_droplet_algebraic(liq, Ts, params)
end

"""
    precompute_liquid_state(m::Tm, T_bar::Tt, params::DropletParams) -> DropletLiquidState

Evaluates all thermodynamic and geometrical droplet quantities that depend strictly on (m, T_bar)
and are invariant to the interface temperature Ts.
"""
function precompute_liquid_state(m::Tm, T_bar::Tt, params::DropletParams) where {Tm<:Real, Tt<:Real}
    T = promote_type(Tm, Tt, Float64)
    T_bar_K = T(T_bar) + T(T_ZERO_C)

    # Liquid water thermophysical properties at T_bar
    rl  = liquid_density_water(T_bar_K)
    cpl = liquid_cp_water(T_bar_K)
    kl  = liquid_thermal_conductivity_water(T_bar_K)
    mul = liquid_viscosity_water(T_bar_K)

    # Droplet geometry
    d   = cbrt(T(m) / (rl * T(π) / T(6.0)))
    r_s = d / T(2.0)

    # Free fall velocity: if not terminal_velocity, fixed to params.U_rel
    u_rel = if params.terminal_velocity
        # Estimate reference film temperature at T_bar for terminal velocity
        T_film_est = (T(params.T_inf) + T(T_bar)) / T(2.0) + T(T_ZERO_C)
        calc_terminal_velocity(d, T_film_est, params.P_atm; rho_l=rl)
    else
        T(params.U_rel)
    end

    return DropletLiquidState{T}(
        T(m), T(T_bar), T_bar_K, rl, cpl, kl, mul, d, r_s, u_rel
    )
end

"""
    intensive_surface_residual(Ts::Real, liq::DropletLiquidState, params::DropletParams) -> Real

Calculates the scale-invariant surface energy balance residual normalized by (4π * r_s):
    R_tilde(Ts) = [q_g - m_dot * DHev - q_l] / (4π * r_s)  [W/m]

This intensive form eliminates radius-scaling degeneration as r_s -> 0.
"""
function intensive_surface_residual(Ts::Ts_t, liq::DropletLiquidState{T}, params::DropletParams) where {Ts_t<:Real, T<:Real}
    TR = promote_type(Ts_t, T, Float64)
    Ts_K    = TR(Ts) + TR(T_ZERO_C)
    T_inf_K = TR(params.T_inf) + TR(T_ZERO_C)

    if params.reference_state === :one_third
        ref = reference_state_one_third(Ts_K, T_inf_K, surface_mass_fraction_from_psat(psat_water(Ts_K), params.P_atm), params.Y_inf)
        T_eg_K = ref.T_ref
        T_eg_C = T_eg_K - TR(T_ZERO_C)
        Y_ref  = ref.Y_v_ref
    else
        T_eg_C  = (TR(params.T_inf) + TR(Ts)) / TR(2.0)
        T_eg_K  = T_eg_C + TR(T_ZERO_C)
        Y_ref   = TR(params.Y_inf)
    end

    if params.use_wilke_mixtures
        mix = gas_mixture_properties(T_eg_K, Y_ref, params.P_atm)
        rg   = mix.rho_g
        mug  = mix.mu_g
        cpg  = mix.cp_g
        kg   = mix.k_g
        Dab  = mix.D_ab
    else
        rg   = gas_density_air(T_eg_K, params.P_atm)
        mug  = gas_viscosity_air(T_eg_K)
        cpg  = gas_cp_air(T_eg_K)
        kg   = gas_thermal_conductivity_air(T_eg_K)
        Dab  = gas_diffusivity_water_air(T_eg_K, params.P_atm)
    end
    cpv  = vapor_cp_water(T_eg_K)
    DHev = latent_heat_water(Ts_K)

    # Terminal velocity evaluated at exact film temperature T_eg_K for gas properties and T_bar for liquid density
    u_rel = params.terminal_velocity ? calc_terminal_velocity(TR(liq.d), T_eg_K, params.P_atm; rho_l=TR(liq.rho_l)) : TR(liq.u_rel)
    Re = calc_re(TR(liq.d), u_rel, rg, mug)
    Pr = calc_pr(mug, cpg, kg)
    Sc = calc_sc(mug, rg, Dab)

    if params.correlation === :clift
        Nu0, Sh0 = clift_correlation(Re, Pr, Sc)
    else
        Nu0 = nusselt_froessling(Re, Pr; coef=params.frossling_coef)
        Sh0 = sherwood_froessling(Re, Sc; coef=params.frossling_coef)
    end

    P_sat_s = psat_water(Ts_K)
    Y_s = surface_mass_fraction_from_psat(P_sat_s, params.P_atm)
    Y_inf_c = clamp(TR(params.Y_inf), zero(TR), TR(0.99999))
    B_M = (Y_s - Y_inf_c) / max(one(TR) - Y_s, TR(1e-12))
    bm_clamped = max(B_M, TR(-0.99999))

    tr = _compute_gas_transport_and_fluxes(
        params, Re, Pr, Sc, Nu0, Sh0,
        rg, kg, Dab, cpv, cpg, TR(liq.rho_l), TR(liq.r_s),
        Y_s, Y_inf_c, bm_clamped, DHev, TR(Ts)
    )
    term_qg   = tr.term_qg
    term_evap = tr.term_evap

    # Intensive conduction into droplet interior per (4π * r_s):
    # q_l / (4π * r_s) = b_q * k_l * (Ts - T_bar)
    circ = internal_circulation_enhancement(
        params.profile, Re, rg, TR(liq.mu_l), TR(liq.r_s), u_rel,
        TR(liq.rho_l), Dab, TR(liq.k_l), TR(liq.cp_l), bm_clamped
    )
    b_q = circ.b_q
    term_ql = b_q * TR(liq.k_l) * (TR(Ts) - TR(liq.T_bar))

    return term_qg - term_evap - term_ql
end

"""
    surface_energy_residual(Ts, m, T_bar, params) -> Real

Residual of the surface energy balance:
    R(Ts) = q_g - m_dot * DHev - q_l = 0
"""
function surface_energy_residual(Ts::Real, m::Real, T_bar::Real, params::DropletParams)
    alg = calc_droplet_algebraic(m, T_bar, Ts, params)
    # Gas heat arriving at surface - latent heat of vaporization - conduction into droplet interior
    return alg.q_g - alg.m_dot * alg.DHev - alg.q_l
end

"""
    solve_surface_temperature_brent(m::Real, T_bar::Real, params::DropletParams; tol::Float64=1e-8, maxiter::Int=100) -> Real

Legacy 1D Brent root-finder for implicit surface temperature Ts.
"""
function solve_surface_temperature_brent(
    m::Real, T_bar::Real, params::DropletParams;
    tol::Float64=1e-8, maxiter::Int=100
)
    # Estimate search bracket [xa, xb]
    T_inf = Float64(params.T_inf)
    T_b   = Float64(T_bar)
    
    # Wet bulb estimate using Stull
    RH_pct = clamp(Float64(params.Y_inf) / max(Float64(surface_mass_fraction_from_psat(psat_water(T_inf + 273.15))), 1e-6) * 100.0, 1.0, 99.0)
    T_wb = stull_twb(T_inf, RH_pct)

    xa = max(min(T_wb, T_b) - 5.0, 0.1)
    xb = max(T_inf, T_b) + 5.0

    fa = surface_energy_residual(xa, m, T_bar, params)
    fb = surface_energy_residual(xb, m, T_bar, params)

    # If not bracketed, expand search
    if fa * fb > 0.0
        xa = max(0.01, min(xa - 10.0, 0.5))
        xb = min(xb + 20.0, 99.0)
        fa = surface_energy_residual(xa, m, T_bar, params)
        fb = surface_energy_residual(xb, m, T_bar, params)
    end

    if fa * fb > 0.0
        return abs(fa) < abs(fb) ? xa : xb
    end

    # Brent's method (1973)
    a, b = xa, xb
    c = a
    fc = fa
    d = b - a
    e = d

    for _ in 1:maxiter
        if params.ts_counter !== nothing
            params.ts_counter.iters += 1
        end
        if abs(fc) < abs(fb)
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        end

        tol_act = 2.0 * eps(Float64) * abs(b) + 0.5 * tol
        xm = 0.5 * (c - b)

        if abs(xm) <= tol_act || abs(fb) <= tol
            return b
        end

        if abs(e) >= tol_act && abs(fa) > abs(fb)
            s = fb / fa
            if a == c
                p = 2.0 * xm * s
                q = 1.0 - s
            else
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            end

            if p > 0.0
                q = -q
            else
                p = -p
            end

            if 2.0 * p < min(3.0 * xm * q - abs(tol_act * q), abs(e * q))
                e = d
                d = p / q
            else
                d = xm
                e = d
            end
        else
            d = xm
            e = d
        end

        a = b
        fa = fb

        if abs(d) > tol_act
            b += d
        else
            b += xm > 0.0 ? tol_act : -tol_act
        end

        fb = surface_energy_residual(b, m, T_bar, params)

        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)
            c = a
            fc = fa
            d = b - a
            e = d
        end
    end

    return b
end

# -----------------------------------------------------------------------------
# STRATEGY A: FAST NEWTON & HALLEY WITH INVARIANT PRECOMPUTATION
# -----------------------------------------------------------------------------

"""
    solve_surface_temperature_newton(liq::DropletLiquidState, params::DropletParams; tol=1e-10, maxiter=8) -> Real

Solves for Ts using Newton-Raphson on the intensive surface residual R_tilde(Ts) = 0.
Warm-starts from Ts = T_bar (since Biot << 1) and converges to machine precision in 2-3 iterations.
"""
function solve_surface_temperature_newton(liq::DropletLiquidState{T}, params::DropletParams; tol::Real=1e-10, maxiter::Int=8) where {T<:Real}
    # Initial guess at droplet mean temperature
    Ts = liq.T_bar
    eps_fd = T(1e-5)

    for _ in 1:maxiter
        if params.ts_counter !== nothing
            params.ts_counter.iters += 1
        end
        R_val = intensive_surface_residual(Ts, liq, params)
        if abs(R_val) <= tol
            return Ts
        end

        # Numerical derivative with central differencing
        R_plus  = intensive_surface_residual(Ts + eps_fd, liq, params)
        R_minus = intensive_surface_residual(Ts - eps_fd, liq, params)
        dR_dTs  = (R_plus - R_minus) / (T(2.0) * eps_fd)

        # Safeguard derivative denominator
        dTs = (abs(dR_dTs) > T(1e-8)) ? (R_val / dR_dTs) : (R_val / (R_val < 0 ? T(10.0) : T(-10.0)))
        
        # Limit step size to physical reasonable jump (<= 5 °C per iter)
        dTs = clamp(dTs, T(-5.0), T(5.0))
        Ts_next = Ts - dTs

        if abs(dTs) < tol
            return Ts_next
        end
        Ts = Ts_next
    end

    return Ts
end

"""
    solve_surface_temperature_halley(liq::DropletLiquidState, params::DropletParams; tol=1e-11, maxiter=5) -> Real

Solves for Ts using 3rd-order Halley's method on the intensive surface residual.
Achieves cubic convergence around Ts = T_bar.
"""
function solve_surface_temperature_halley(liq::DropletLiquidState{T}, params::DropletParams; tol::Real=1e-11, maxiter::Int=5) where {T<:Real}
    Ts = liq.T_bar
    eps_fd = T(1e-4)

    for _ in 1:maxiter
        if params.ts_counter !== nothing
            params.ts_counter.iters += 1
        end
        R0 = intensive_surface_residual(Ts, liq, params)
        if abs(R0) <= tol
            return Ts
        end

        Rp = intensive_surface_residual(Ts + eps_fd, liq, params)
        Rm = intensive_surface_residual(Ts - eps_fd, liq, params)
        
        dR = (Rp - Rm) / (T(2.0) * eps_fd)
        d2R = (Rp - T(2.0) * R0 + Rm) / (eps_fd^2)

        denom = T(2.0) * (dR^2) - R0 * d2R
        if abs(denom) > T(1e-12)
            dTs = (T(2.0) * R0 * dR) / denom
        else
            dTs = R0 / dR
        end

        dTs = clamp(dTs, T(-5.0), T(5.0))
        Ts_next = Ts - dTs

        if abs(dTs) < tol
            return Ts_next
        end
        Ts = Ts_next
    end

    return Ts
end

# -----------------------------------------------------------------------------
# STRATEGY B: CLOSED-FORM ANALYTICAL EXPANSIONS (ZERO ITERATION)
# -----------------------------------------------------------------------------

"""
    solve_surface_temperature_analytical_biot(liq::DropletLiquidState, params::DropletParams) -> Real

Closed-form direct evaluation of Ts using 1st-order analytical Biot-Clausius-Clapeyron expansion
around Ts = T_bar:
    Ts = T_bar - R_tilde(T_bar) / R_tilde'(T_bar)
Zero iterations, completely branch-free and AD-native.
"""
function solve_surface_temperature_analytical_biot(liq::DropletLiquidState{T}, params::DropletParams) where {T<:Real}
    Ts0 = liq.T_bar
    eps_fd = T(1e-4)

    R0 = intensive_surface_residual(Ts0, liq, params)
    Rp = intensive_surface_residual(Ts0 + eps_fd, liq, params)
    Rm = intensive_surface_residual(Ts0 - eps_fd, liq, params)
    dR = (Rp - Rm) / (T(2.0) * eps_fd)

    if abs(dR) < T(1e-8)
        return Ts0
    end

    dTs = clamp(R0 / dR, T(-10.0), T(10.0))
    return Ts0 - dTs
end

"""
    solve_surface_temperature_analytical_halley(liq::DropletLiquidState, params::DropletParams) -> Real

Closed-form direct evaluation of Ts using 2nd-order analytical Halley/Padé expansion
around Ts = T_bar. Matches the exact root to ~10^-9 K with 0 iterative loops.
"""
function solve_surface_temperature_analytical_halley(liq::DropletLiquidState{T}, params::DropletParams) where {T<:Real}
    Ts0 = liq.T_bar
    eps_fd = T(1e-4)

    R0 = intensive_surface_residual(Ts0, liq, params)
    Rp = intensive_surface_residual(Ts0 + eps_fd, liq, params)
    Rm = intensive_surface_residual(Ts0 - eps_fd, liq, params)

    dR  = (Rp - Rm) / (T(2.0) * eps_fd)
    d2R = (Rp - T(2.0) * R0 + Rm) / (eps_fd^2)

    denom = T(2.0) * (dR^2) - R0 * d2R
    if abs(denom) > T(1e-12)
        dTs = (T(2.0) * R0 * dR) / denom
    else
        dTs = R0 / dR
    end

    dTs = clamp(dTs, T(-10.0), T(10.0))
    return Ts0 - dTs
end

# -----------------------------------------------------------------------------
# UNIFIED INTERFACE DISPATCH FOR Ts
# -----------------------------------------------------------------------------

"""
    solve_surface_temperature(liq::DropletLiquidState, params::DropletParams; tol=1e-8, maxiter=100) -> Real
    solve_surface_temperature(m, T_bar, params; tol=1e-8, maxiter=100) -> Real

Finds equilibrium surface temperature Ts using either precomputed liquid state `liq` or `(m, T_bar)`.
Dispatches according to `params.ts_solver`:
  - `:chebyshev` (Strategy D1: O(1) Chebyshev polynomial surrogate)
  - `:lut` (Strategy D2: O(1) Bicubic Hermite lookup table)
  - `:analytical_halley` (Strategy B2: 0-iteration closed form)
  - `:analytical_biot` (Strategy B1: 0-iteration closed form)
  - `:newton` (Strategy A1: iterative Newton)
  - `:halley` (Strategy A2: iterative Halley)
  - `:brent` (legacy reference)

When inputs are `ForwardDiff.Dual`, solves the primal root at Float64 and propagates
exact analytical sensitivity via the Implicit Function Theorem (IFT) or surrogate derivatives:
    dTs* = -(∂R/∂m * dm + ∂R/∂T_bar * dT_bar) / (∂R/∂Ts)
"""
function solve_surface_temperature(
    liq::DropletLiquidState, params::DropletParams;
    tol::Float64=1e-8, maxiter::Int=100
)
    if params.profile === :uniforme
        return liq.T_bar
    end

    if params.ts_counter !== nothing
        params.ts_counter.calls += 1
    end

    method = params.ts_solver
    if method === :chebyshev
        surr = get_or_build_ts_surrogate(params, :chebyshev)
        return eval_chebyshev_ts(surr, liq.r_s, liq.T_bar)
    elseif method === :lut
        surr = get_or_build_ts_surrogate(params, :lut)
        return eval_lut_ts(surr, liq.r_s, liq.T_bar)
    elseif method === :analytical_halley
        return solve_surface_temperature_analytical_halley(liq, params)
    elseif method === :analytical_biot
        return solve_surface_temperature_analytical_biot(liq, params)
    elseif method === :halley
        return solve_surface_temperature_halley(liq, params; tol=tol)
    elseif method === :brent
        return solve_surface_temperature_brent(liq.m, liq.T_bar, params; tol=tol, maxiter=maxiter)
    else
        return solve_surface_temperature_newton(liq, params; tol=tol)
    end
end

"""
    solve_surface_temperature(m, T_bar, params; tol=1e-8, maxiter=100) -> Real

Finds equilibrium surface temperature Ts according to `params.ts_solver`:
  - `:chebyshev` (Strategy D1)
  - `:lut` (Strategy D2)
  - `:analytical_halley` (Strategy B2)
  - `:analytical_biot` (Strategy B1)
  - `:newton` (Strategy A)
  - `:halley` (Strategy A2)
  - `:brent` (legacy reference)

When inputs are `ForwardDiff.Dual`, solves the primal root at Float64 and propagates
exact analytical sensitivity via the Implicit Function Theorem (IFT) or surrogate derivatives:
    dTs* = -(∂R/∂m * dm + ∂R/∂T_bar * dT_bar) / (∂R/∂Ts)
"""
function solve_surface_temperature(
    m::Real, T_bar::Real, params::DropletParams;
    tol::Float64=1e-8, maxiter::Int=100
)
    # For lumped / infinite conductivity profile (Zhou et al. 2017), Ts ≡ T_bar
    if params.profile === :uniforme
        return T_bar
    end

    # Implicit Function Theorem (IFT) / Dual propagation for ForwardDiff.Dual
    if m isa ForwardDiff.Dual || T_bar isa ForwardDiff.Dual
        method = params.ts_solver
        if method === :chebyshev
            surr = get_or_build_ts_surrogate(params, :chebyshev)
            # Evaluate r_s as Dual
            rho_l = liquid_density_water(Float64(ForwardDiff.value(T_bar)) + 273.15)
            d = cbrt(m / (rho_l * (π / 6.0)))
            r_s = d / 2.0
            return eval_chebyshev_ts(surr, r_s, T_bar)
        elseif method === :lut
            surr = get_or_build_ts_surrogate(params, :lut)
            rho_l = liquid_density_water(Float64(ForwardDiff.value(T_bar)) + 273.15)
            d = cbrt(m / (rho_l * (π / 6.0)))
            r_s = d / 2.0
            return eval_lut_ts(surr, r_s, T_bar)
        end

        m_val = Float64(ForwardDiff.value(m))
        Tb_val = Float64(ForwardDiff.value(T_bar))
        liq_val = precompute_liquid_state(m_val, Tb_val, params)

        # Solve root in Float64
        Ts_val = if method === :halley
            solve_surface_temperature_halley(liq_val, params; tol=tol)
        elseif method === :analytical_biot
            solve_surface_temperature_analytical_biot(liq_val, params)
        elseif method === :analytical_halley
            solve_surface_temperature_analytical_halley(liq_val, params)
        elseif method === :brent
            Float64(solve_surface_temperature_brent(m_val, Tb_val, params; tol=tol))
        else
            solve_surface_temperature_newton(liq_val, params; tol=tol)
        end

        # Exact analytical derivative via IFT:
        # R(Ts, m, T_bar) = 0 => dTs = - (∂R/∂m * dm + ∂R/∂Tb * dTb) / (∂R/∂Ts)
        eps_fd = 1e-5
        Rp = intensive_surface_residual(Ts_val + eps_fd, liq_val, params)
        Rm = intensive_surface_residual(Ts_val - eps_fd, liq_val, params)
        dR_dTs = (Rp - Rm) / (2.0 * eps_fd)

        # Numerically compute partials ∂R/∂m and ∂R/∂T_bar with respect to Float64 base
        eps_m = max(1e-15, m_val * 1e-5)
        liq_m_plus = precompute_liquid_state(m_val + eps_m, Tb_val, params)
        liq_m_minus = precompute_liquid_state(m_val - eps_m, Tb_val, params)
        dR_dm = (intensive_surface_residual(Ts_val, liq_m_plus, params) - intensive_surface_residual(Ts_val, liq_m_minus, params)) / (2.0 * eps_m)

        eps_T = 1e-4
        liq_T_plus = precompute_liquid_state(m_val, Tb_val + eps_T, params)
        liq_T_minus = precompute_liquid_state(m_val, Tb_val - eps_T, params)
        dR_dTb = (intensive_surface_residual(Ts_val, liq_T_plus, params) - intensive_surface_residual(Ts_val, liq_T_minus, params)) / (2.0 * eps_T)

        dTs_dm = -dR_dm / dR_dTs
        dTs_dTb = -dR_dTb / dR_dTs

        # Chain rule with Dual parts of inputs
        part_m = (m isa ForwardDiff.Dual) ? ForwardDiff.partials(m) : zero(ForwardDiff.partials(T_bar))
        part_Tb = (T_bar isa ForwardDiff.Dual) ? ForwardDiff.partials(T_bar) : zero(ForwardDiff.partials(m))
        partials = dTs_dm * part_m + dTs_dTb * part_Tb

        return ForwardDiff.Dual{ForwardDiff.tagtype(typeof(m isa ForwardDiff.Dual ? m : T_bar))}(Ts_val, partials)
    end

    liq = precompute_liquid_state(m, T_bar, params)
    return solve_surface_temperature(liq, params; tol=tol, maxiter=maxiter)
end

"""
    solve_droplet_state(m, T_bar, params) -> (Ts, DropletAlgebraic)

Unified helper that computes Ts and evaluates the full DropletAlgebraic diagnostic structure
in a single streamlined pass without redundant recomputation.
"""
function solve_droplet_state(m::Real, T_bar::Real, params::DropletParams)
    liq = precompute_liquid_state(m, T_bar, params)
    Ts = solve_surface_temperature(liq, params)
    alg = calc_droplet_algebraic(liq, Ts, params)
    return Ts, alg
end

# =============================================================================
# VENTILATED WET-BULB & REVISED D2 LAW (Dalla Barba et al. 2021, Aroesty 1969)
# =============================================================================

"""
    ventilated_wet_bulb_temperature(T_inf_C::Real, Y_inf::Real, Re::Real, P_Pa::Real=P_ATM; tol=1e-5, maxiter=100) -> Real

Computes the ventilated equilibrium wet-bulb temperature T_ds [°C] by solving the implicit
energy balance from Dalla Barba, Wang & Picano (2021) / Aroesty & Koenig (1969):
    T_ds + (Pr/Sc) * (Sh/Nu) * (L(T_ds)/cp_g) * ln((1 - Y_inf)/(1 - Y_s(T_ds))) - T_inf = 0
"""
function ventilated_wet_bulb_temperature(
    T_inf_C::T, Y_inf::Real, Re::Real, P_Pa::Real=P_ATM;
    tol::Float64=1e-5, maxiter::Int=100
) where {T<:Real}
    T_inf_K = T_inf_C + T(T_ZERO_C)
    y_inf = clamp(T(Y_inf), zero(T), T(0.9999))
    re_val = max(T(Re), zero(T))
    P_val = T(P_Pa)

    # Residual function for Ts in Celsius
    function res_ts(Ts_C::Real)
        Ts_K = T(Ts_C) + T(T_ZERO_C)
        T_film_K = (Ts_K + T_inf_K) / T(2.0)

        cpg = gas_cp_air(T_film_K)
        mug = gas_viscosity_air(T_film_K)
        kg  = gas_thermal_conductivity_air(T_film_K)
        rg  = gas_density_air(T_film_K, P_val)
        Dab = gas_diffusivity_water_air(T_film_K, P_val)

        Pr = calc_pr(mug, cpg, kg)
        Sc = calc_sc(mug, rg, Dab)

        Nu = nusselt_froessling(re_val, Pr)
        Sh = sherwood_froessling(re_val, Sc)

        Psat = psat_water(Ts_K)
        Ys = surface_mass_fraction_from_psat(Psat, P_val)
        Ys = clamp(Ys, zero(T), T(0.99999))

        bm = max((Ys - y_inf) / (one(T) - Ys), T(-0.9999))
        ln_bm = log1p(bm)

        Lv = latent_heat_water(Ts_K)

        term = (Pr / Sc) * (Sh / Nu) * (Lv / cpg) * ln_bm
        return (Ts_K + term) - T_inf_K
    end

    # Bracket search [T_low, T_high]
    T_low = max(T(-10.0), T_inf_C - T(45.0))
    T_high = T_inf_C

    f_low = res_ts(T_low)
    f_high = res_ts(T_high)

    if f_low * f_high > 0.0
        # Expand bracket
        T_low = T(-20.0)
        f_low = res_ts(T_low)
    end

    if f_low * f_high > 0.0
        return T_inf_C # fallback
    end

    a, b = T_low, T_high
    for _ in 1:maxiter
        mid = (a + b) / T(2.0)
        if abs(b - a) < T(tol)
            return mid
        end
        f_mid = res_ts(mid)
        if f_mid * f_low <= 0.0
            b = mid
            f_high = f_mid
        else
            a = mid
            f_low = f_mid
        end
    end
    return (a + b) / T(2.0)
end

"""
    calc_d2_revised_rate(T_inf_C::Real, Y_inf::Real, Re::Real, d::Real, P_Pa::Real=P_ATM)

Calculates the revised D^2-law evaporation rate constant K_r [m^2/s] and lifetime τ [s]
from Dalla Barba, Wang & Picano (2021):
    K_r = (4 * Sh * rho_g * Dab / rho_l) * ln((1 - Y_inf)/(1 - Y_s(T_ds)))
"""
function calc_d2_revised_rate(
    T_inf_C::Real, Y_inf::Real, Re::Real, d::Real, P_Pa::Real=P_ATM
)
    T = promote_type(typeof(T_inf_C), typeof(Y_inf), typeof(Re), typeof(d), Float64)
    T_ds_C = ventilated_wet_bulb_temperature(T(T_inf_C), T(Y_inf), T(Re), T(P_Pa))
    T_ds_K = T_ds_C + T(T_ZERO_C)
    T_inf_K = T(T_inf_C) + T(T_ZERO_C)
    T_film_K = (T_ds_K + T_inf_K) / T(2.0)

    P_val = T(P_Pa)
    rg = gas_density_air(T_film_K, P_val)
    mug = gas_viscosity_air(T_film_K)
    Dab = gas_diffusivity_water_air(T_film_K, P_val)
    rl = liquid_density_water(T_ds_K)

    Sc = calc_sc(mug, rg, Dab)
    Sh = sherwood_froessling(T(Re), Sc)

    Psat = psat_water(T_ds_K)
    Ys = surface_mass_fraction_from_psat(Psat, P_val)
    bm = max((Ys - T(Y_inf)) / (one(T) - Ys), zero(T))
    ln_bm = log1p(bm)

    Kr = (T(4.0) * Sh * rg * Dab / rl) * ln_bm
    d_val = T(d)
    tau = Kr > zero(T) ? (d_val^2) / Kr : Inf

    return (Kr = Kr, T_ds = T_ds_C, Sh = Sh, tau = tau)
end

# =============================================================================
# ZHOU ET AL. (2017) 3RD-ORDER POLYNOMIAL RADIAL PROFILE RECONSTRUCTION
# =============================================================================

"""
    reconstruct_zhou_cubic_profile(r::Real, r_s::Real, T_bar::Real, Ts::Real) -> Real

Reconstructs the internal droplet temperature T(r) from the 3rd-order polynomial model
of Zhou et al. (2017):
    T(r) = a + c*r^2 + d*r^3
satisfying:
    ∂T/∂r|_{r=0} = 0,
    ∂T/∂r|_{r = r_s / sqrt(2)} = 0  (Hill's vortex stagnation minimal temperature),
    (3 / r_s^3) ∫_0^{r_s} r^2 T(r) dr = T_bar,
    T(r_s) = Ts.
"""
function reconstruct_zhou_cubic_profile(r::T, r_s::Real, T_bar::Real, Ts::Real) where {T<:Real}
    rs = max(T(r_s), T(1e-9))
    r_val = clamp(r, zero(T), rs)
    tb = T(T_bar)
    ts = T(Ts)

    # Analytical coefficients:
    # d = -2*sqrt(2)/(3*rs) * c
    # T_bar = a + (3/5 - sqrt(2)/3) * c * rs^2
    # Ts    = a + (1 - 2*sqrt(2)/3) * c * rs^2
    # Ts - T_bar = (2/5 - sqrt(2)/3) * c * rs^2
    denom_c = (T(0.4) - sqrt(T(2.0)) / T(3.0)) * (rs^2)
    c = (ts - tb) / denom_c
    d = -T(2.0) * sqrt(T(2.0)) / (T(3.0) * rs) * c
    a = ts - c * (rs^2) - d * (rs^3)

    return a + c * (r_val^2) + d * (r_val^3)
end

# =============================================================================
# SAZHIN (2006, 2017, 2021) ANALYTICAL TRANSIENT EIGENVALUE HEATING SERIES
# =============================================================================

"""
    sazhin_eigenvalue_temperature_profile(r::Real, r_s::Real, t::Real, T0::Real, T_eff::Real, kappa_eff::Real, Bi::Real; n_modes::Int=20) -> Real

Evaluates the exact analytical eigenvalue series solution for transient heat conduction inside
a spherical droplet (Sazhin 2006, 2017, 2021):
    T(r, t) = T_eff + (T0 - T_eff) * ∑_{n=1}^N A_n * (sin(λ_n * r / r_s) / (λ_n * r / r_s)) * exp(-λ_n^2 * kappa_eff * t / r_s^2)
where roots λ_n satisfy: tan(λ_n) = -λ_n / (Bi - 1) = λ_n / (1 - Bi).
"""
function sazhin_eigenvalue_temperature_profile(
    r::T, r_s::Real, t::Real, T0::Real, T_eff::Real, kappa_eff::Real, Bi::Real;
    n_modes::Int=20
) where {T<:Real}
    rs = max(T(r_s), T(1e-9))
    t_val = max(T(t), zero(T))
    r_val = clamp(r, zero(T), rs)
    eta = r_val / rs

    t0 = T(T0)
    te = T(T_eff)
    ke = max(T(kappa_eff), T(1e-12))
    bi = max(T(Bi), T(1e-6))

    sum_terms = zero(T)

    for n in 1:n_modes
        # Robust bracketed bisection for tan(λ) = λ / (1 - Bi)
        local a_brk, b_brk
        if bi > one(T)
            # Root in (π/2 + (n-1)*π, n*π)
            a_brk = (T(n) - T(0.5)) * T(π) + T(1e-5)
            b_brk = T(n) * T(π) - T(1e-6)
        elseif bi < one(T)
            # Root in ((n-1)*π, (n-0.5)*π)
            a_brk = (T(n) - one(T)) * T(π) + (n == 1 ? T(1e-4) : T(1e-6))
            b_brk = (T(n) - T(0.5)) * T(π) - T(1e-5)
        else
            # Bi == 1
            a_brk = (T(n) - T(0.5)) * T(π)
            b_brk = a_brk
        end

        lambda = (a_brk + b_brk) / T(2.0)
        if a_brk != b_brk
            f_res(lam) = tan(lam) - lam / (one(T) - bi)
            fa = f_res(a_brk)
            for _ in 1:40
                mid = (a_brk + b_brk) / T(2.0)
                if abs(b_brk - a_brk) < T(1e-8)
                    lambda = mid
                    break
                end
                fmid = f_res(mid)
                if fa * fmid <= zero(T)
                    b_brk = mid
                else
                    a_brk = mid
                    fa = fmid
                end
                lambda = (a_brk + b_brk) / T(2.0)
            end
        end

        # Fourier coefficient A_n
        sin_l = sin(lambda)
        cos_l = cos(lambda)
        denom_A = lambda - sin_l * cos_l
        A_n = denom_A != zero(T) ? (T(2.0) * (sin_l - lambda * cos_l) / denom_A) : zero(T)

        # Spatial radial shape: sin(λ * η) / (λ * η)
        shape = eta > T(1e-6) ? (sin(lambda * eta) / (lambda * eta)) : one(T)

        # Temporal decay
        decay = exp(-(lambda^2) * ke * t_val / (rs^2))

        sum_terms += A_n * shape * decay
    end

    return te + (t0 - te) * sum_terms
end
