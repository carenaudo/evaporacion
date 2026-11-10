"""
    Kinetic droplet evaporation solvers and Knudsen-layer non-equilibrium formulations.
    Based on Vlasov (2021), Zhao & Nadal (2023), and Akkuş (2020).
"""

using ..Common: R_UNIVERSAL, P_ATM, psat_water, gas_diffusivity_water_air, liquid_density_water
using OrdinaryDiffEq

"""
    kelvin_vapor_pressure_ratio(r, T_s; rho_l=997.0, surface_tension=0.0728, molar_mass=0.018015)

Compute the Kelvin curvature enhancement factor for equilibrium vapor pressure:
```math
\\frac{p_{\\text{sat,curv}}}{p_{\\text{sat}}} = \\exp\\left( \\frac{2 \\sigma M_v}{\\rho_l R_g T_s r} \\right)
```
Valid for sub-micron and nanoscale droplets where surface curvature raises the chemical potential.
"""
function kelvin_vapor_pressure_ratio(r::Real, T_s::Real;
                                     rho_l::Real=997.0,
                                     surface_tension::Real=0.0728,
                                     molar_mass::Real=0.018015)
    r_val = float(r)
    if r_val <= 0.0
        return 1.0
    end
    arg = (2.0 * surface_tension * molar_mass) / (rho_l * R_UNIVERSAL * T_s * r_val)
    # Cap argument to prevent numerical overflow for near-zero radii
    return exp(clamp(arg, 0.0, 50.0))
end

"""
    psat_kelvin(r, T_s; rho_l=997.0, surface_tension=0.0728, molar_mass=0.018015)

Calculate Kelvin-corrected saturation vapor pressure [Pa] at droplet radius `r` and temperature `T_s`.
"""
function psat_kelvin(r::Real, T_s::Real;
                     rho_l::Real=997.0,
                     surface_tension::Real=0.0728,
                     molar_mass::Real=0.018015)
    p_flat = psat_water(T_s)
    k_factor = kelvin_vapor_pressure_ratio(r, T_s; rho_l=rho_l, surface_tension=surface_tension, molar_mass=molar_mass)
    return p_flat * k_factor
end

"""
    mean_thermal_speed(T; molar_mass=0.018015)

Calculate mean Maxwellian molecular thermal speed [m/s]:
```math
\\upsilon = \\sqrt{\\frac{8 R_g T}{\\pi M_v}}
```
"""
function mean_thermal_speed(T::Real; molar_mass::Real=0.018015)
    return sqrt((8.0 * R_UNIVERSAL * float(T)) / (pi * float(molar_mass)))
end

"""
    vlasov_kinetic_length(D, alpha, v)

Calculate the characteristic kinetic transition length [m] from Vlasov (2021):
```math
l_k = \\frac{4 D}{\\alpha \\upsilon}
```
When droplet radius `r >> l_k`, evaporation is diffusion-controlled (classical d^2 law).
When `r << l_k`, evaporation is interface kinetics-controlled (linear d law).
"""
function vlasov_kinetic_length(D::Real, alpha::Real, v::Real)
    return (4.0 * float(D)) / (float(alpha) * float(v))
end

"""
    vlasov_evaporation_rate(r, D, delta_c, alpha, v)

Compute the diffusion-kinetic molar evaporation rate [mol/s] according to Vlasov (2021) Eq. (18):
```math
I = \\frac{4 \\pi r D \\Delta c}{1 + \\frac{4 D}{\\alpha \\upsilon r}} = \\frac{4 \\pi r D \\Delta c}{1 + l_k / r}
```
"""
function vlasov_evaporation_rate(r::Real, D::Real, delta_c::Real, alpha::Real, v::Real)
    r_val = float(r)
    if r_val <= 0.0 || delta_c <= 0.0
        return 0.0
    end
    lk = vlasov_kinetic_length(D, alpha, v)
    return (4.0 * pi * r_val * D * delta_c) / (1.0 + lk / r_val)
end

"""
    vlasov_radius_trajectory(t, r0, D, delta_c, rho_l, molar_mass, alpha, v)

Exact analytical droplet radius regression [m] at time `t` [s] from Vlasov (2021) Eq. (20):
```math
R(t) = -l_k + \\sqrt{(r_0 + l_k)^2 - \\frac{2 D M_v}{\\rho_l} \\Delta c\\, t}
```
"""
function vlasov_radius_trajectory(t::Real, r0::Real, D::Real, delta_c::Real,
                                  rho_l::Real, molar_mass::Real, alpha::Real, v::Real)
    lk = vlasov_kinetic_length(D, alpha, v)
    term_diff = (2.0 * D * molar_mass * delta_c * float(t)) / float(rho_l)
    radicand = (float(r0) + lk)^2 - term_diff
    if radicand <= lk^2
        return 0.0
    end
    return -lk + sqrt(radicand)
end

"""
    vlasov_lifetime(r0, D, delta_c, rho_l, molar_mass, alpha, v)

Exact analytical droplet lifetime [s] from Vlasov (2021), extending the Maxwellian d^2 lifetime:
```math
\\tau_{\\text{Vlasov}} = \\frac{\\rho_l (r_0^2 + 2 r_0 l_k)}{2 D M_v \\Delta c} = \\tau_{\\text{Maxwell}} \\left(1 + \\frac{2 l_k}{r_0}\\right)
```
"""
function vlasov_lifetime(r0::Real, D::Real, delta_c::Real,
                         rho_l::Real, molar_mass::Real, alpha::Real, v::Real)
    if delta_c <= 0.0 || r0 <= 0.0
        return Inf
    end
    lk = vlasov_kinetic_length(D, alpha, v)
    tau_maxwell = (float(rho_l) * float(r0)^2) / (2.0 * float(D) * float(molar_mass) * float(delta_c))
    return tau_maxwell * (1.0 + (2.0 * lk) / float(r0))
end

"""
    kinetic_peclet_number(r, T_inf, D; molar_mass=0.018015)

Calculate the kinetic Péclet number from Zhao & Nadal (2023):
```math
\\text{Pe}_k = \\frac{u_k r}{\\mathcal{D}_{AB}}, \\qquad u_k = \\sqrt{\\frac{R_g T_\\infty}{2 \\pi M_v}}
```
Measures the ratio between kinetic molecular effusion speed and macroscopic diffusion speed.
"""
function kinetic_peclet_number(r::Real, T_inf::Real, D::Real; molar_mass::Real=0.018015)
    u_k = sqrt((R_UNIVERSAL * float(T_inf)) / (2.0 * pi * float(molar_mass)))
    return (u_k * float(r)) / float(D)
end

"""
    zhao_nadal_lifetime(tau_d2, Pe0; c3=nothing)

Calculates the non-equilibrium extended droplet lifetime from Zhao & Nadal (2023) Eq. (3.8):
```math
\\tau_{\\text{evap}} = \\tau_{D^2} (1 + c_3)
```
where ``c_3 \\propto 1 / \\text{Pe}_0``. If `c3` is not provided, defaults to ``c_3 = 2 / \\text{Pe}_0``.
"""
function zhao_nadal_lifetime(tau_d2::Real, Pe0::Real; c3::Union{Nothing, Real}=nothing)
    c3_val = c3 === nothing ? (2.0 / float(Pe0)) : float(c3)
    return float(tau_d2) * (1.0 + max(0.0, c3_val))
end

"""
    stefan_flow_enhancement(Y_vs, Y_v_inf=0.0)

Stefan flow interfacial blowing correction factor (Akkuş, 2020):
```math
\\Theta_{\\text{Stefan}} = \\frac{1}{1 - Y_{v,s}}
```
Enhances the mass flux when vapor concentration near the surface is non-negligible.
"""
function stefan_flow_enhancement(Y_vs::Real, Y_v_inf::Real=0.0)
    y_s = clamp(float(Y_vs), 0.0, 0.99)
    return 1.0 / (1.0 - y_s)
end

"""
    solve_kinetic_droplet(cfg::KineticDropletConfig{T}; t_span=nothing, r_min=1e-9)

Solve transient sub-micron or microscale droplet evaporation governed by non-equilibrium
Knudsen-layer kinetics, Kelvin curvature shifts, and Stefan blowing.
"""
function solve_kinetic_droplet(cfg::KineticDropletConfig{T};
                                t_span=nothing,
                                r_min::Real=1e-9) where {T<:Real}
    D = gas_diffusivity_water_air(cfg.T_inf, cfg.P_tot)
    v_th = mean_thermal_speed(cfg.T_s; molar_mass=cfg.molar_mass)
    lk = vlasov_kinetic_length(D, cfg.alpha, v_th)

    # Ambient vapor concentration
    p_sat_inf = psat_water(cfg.T_inf)
    p_v_inf = cfg.RH * p_sat_inf
    c_inf = p_v_inf / (R_UNIVERSAL * cfg.T_inf)

    # Base surface vapor concentration (flat interface)
    p_sat_s_flat = psat_water(cfg.T_s)
    c_eq_flat = p_sat_s_flat / (R_UNIVERSAL * cfg.T_s)
    delta_c_flat = max(0.0, c_eq_flat - c_inf)

    # Continuum Maxwellian lifetime for comparison
    tau_cont = delta_c_flat > 0.0 ?
        (cfg.rho_l * cfg.r0^2) / (2.0 * D * cfg.molar_mass * delta_c_flat) : Inf

    tau_est = delta_c_flat > 0.0 ?
        vlasov_lifetime(cfg.r0, D, delta_c_flat, cfg.rho_l, cfg.molar_mass, cfg.alpha, v_th) : 100.0

    t_end = t_span === nothing ? (isfinite(tau_est) ? 1.5 * tau_est : 100.0) : t_span[2]

    # ODE: dr/dt = - (M_v / rho_l) * (c_eq(r) - c_inf) / (4 / (alpha * v) + r / D) * Stefan
    function dr_dt(r_val, p, t)
        if r_val <= r_min
            return 0.0
        end
        # Kelvin effect
        c_eq = if cfg.use_kelvin
            k_ratio = kelvin_vapor_pressure_ratio(r_val, cfg.T_s;
                                                  rho_l=cfg.rho_l,
                                                  surface_tension=cfg.surface_tension,
                                                  molar_mass=cfg.molar_mass)
            (p_sat_s_flat * k_ratio) / (R_UNIVERSAL * cfg.T_s)
        else
            c_eq_flat
        end

        delta_c = c_eq - c_inf
        if delta_c <= 0.0
            return 0.0
        end

        # Stefan flow correction
        stefan_factor = if cfg.use_stefan
            p_s = cfg.use_kelvin ? p_sat_s_flat * kelvin_vapor_pressure_ratio(r_val, cfg.T_s; rho_l=cfg.rho_l, surface_tension=cfg.surface_tension, molar_mass=cfg.molar_mass) : p_sat_s_flat
            Y_vs = (p_s * cfg.molar_mass) / (p_s * cfg.molar_mass + (cfg.P_tot - p_s) * cfg.carrier_molar_mass)
            stefan_flow_enhancement(Y_vs)
        else
            1.0
        end

        denom = (4.0 / (cfg.alpha * v_th)) + (r_val / D)
        rate = -(cfg.molar_mass / cfg.rho_l) * (delta_c / denom) * stefan_factor
        return rate
    end

    condition(u, t, integrator) = u - r_min
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    prob = ODEProblem(dr_dt, float(cfg.r0), (0.0, float(t_end)))
    sol = solve(prob, Tsit5(); callback=cb, reltol=1e-6, abstol=1e-10)

    t_vec = sol.t
    r_vec = clamp.(sol.u, 0.0, float(cfg.r0))
    mass_vec = [(4.0 / 3.0) * pi * r^3 * cfg.rho_l for r in r_vec]
    final_time = t_vec[end]

    pe0 = kinetic_peclet_number(cfg.r0, cfg.T_inf, D; molar_mass=cfg.molar_mass)

    return KineticSolution{T}(
        t_vec,
        r_vec,
        mass_vec,
        final_time,
        tau_cont,
        final_time / tau_cont
    )
end
