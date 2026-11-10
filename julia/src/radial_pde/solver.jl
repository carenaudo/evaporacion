"""
    solver.jl

1D Moving-boundary radial conduction PDE solver for evaporating droplets.
Solves heat conduction on an immobilized Landau coordinate ξ = r / R(t) ∈ [0, 1]
with Marangoni interfacial enhancement, natural convection, and surfactant transport.
"""

"""
    RadialPDEConfig

Configuration for 1D moving-boundary radial conduction PDE.
"""
Base.@kwdef struct RadialPDEConfig
    n_nodes::Int = 31                       # Number of radial grid points
    d0::Float64 = 1.0e-3                    # Initial diameter [m]
    T_drop_0::Float64 = 20.0                # Initial droplet temperature [°C]
    T_air::Float64 = 25.0                   # Ambient air temperature [°C]
    RH::Float64 = 0.5                       # Relative humidity [0, 1]
    with_surfactants::Bool = false          # Enable surfactant concentration & adsorption
    with_marangoni::Bool = true             # Enable Marangoni interfacial flow enhancement
    with_natural_convection::Bool = true    # Enable Rayleigh/Grashof buoyancy enhancement
    m_s_init::Float64 = 1.0e-6              # Solute initial mass [kg]
    m_w_init::Float64 = 1.0e-3              # Water initial mass [kg]
    Gamma_inf::Float64 = 3.0e-6             # Maximum surface excess [mol/m^2]
    k_ads::Float64 = 3.0                    # Adsorption rate constant [m^3/(mol*s)]
    k_des::Float64 = 0.05                   # Desorption rate constant [1/s]
    t_end::Float64 = 1000.0                 # Integration time horizon [s]
end

struct RadialPDESolution
    t::Vector{Float64}
    R::Vector{Float64}
    T_profile::Matrix{Float64}
    m_w::Vector{Float64}
    Gamma::Vector{Float64}
end

"""
    sherwood_natural_convection(R, Ts_K, T_air_K, rho_g, mu_g, Dv) -> Real

Sherwood number with buoyancy / natural convection (Grashof / Rayleigh number).
"""
function sherwood_natural_convection(R::Real, Ts_K::Real, T_air_K::Real, rho_g::Real, mu_g::Real, Dv::Real)
    d = 2.0 * R
    Sc = calc_sc(mu_g, rho_g, Dv)
    T_film = 0.5 * (T_air_K + Ts_K)
    beta = 1.0 / T_film
    deltaT = abs(T_air_K - Ts_K)
    nu_g = mu_g / rho_g
    Gr = GRAVITY * beta * deltaT * (d^3) / (nu_g^2)
    Ra = Gr * Sc
    Sh_Ra = 0.54 * (max(Ra, 0.0)^0.25)
    return (2.0^3 + Sh_Ra^3)^(1.0 / 3.0)
end

"""
    marangoni_factor(R, Ts_K, T_air_K, mu_l, alpha_l) -> Real

Interfacial Marangoni convection evaporation enhancement factor.
"""
function marangoni_factor(R::Real, Ts_K::Real, T_air_K::Real, mu_l::Real, alpha_l::Real)
    gamma_T = 1.5e-4 # |dσ/dT| [N/(m*K)]
    deltaT = abs(T_air_K - Ts_K)
    Ma = gamma_T * R * deltaT / (max(mu_l, 1e-12) * max(alpha_l, 1e-16))
    return 1.0 + 0.2 * (max(Ma, 0.0)^(1.0 / 3.0))
end

"""
    radial_pde_rhs!(dy, y, p::RadialPDEConfig, t)

Semi-discrete Method of Lines for the 1D immobilized Landau coordinate moving-boundary PDE.
State vector:
    y[1] = R(t)
    y[2 : N+1] = T(ξ_i, t) in °C
    y[N+2] = m_w (if with_surfactants)
    y[N+3] = Gamma (if with_surfactants)
"""
function radial_pde_rhs!(dy, y, p::RadialPDEConfig, t)
    N = p.n_nodes
    R = y[1]

    if R <= 1e-7
        dy .= 0.0
        return nothing
    end

    T_nodes = @view y[2:N+1]
    Ts_C = T_nodes[end]
    Ts_K = Ts_C + T_ZERO_C
    T_air_K = p.T_air + T_ZERO_C
    T_film_K = 0.5 * (Ts_K + T_air_K)

    # Physical properties
    rho_g = gas_density_air(T_film_K)
    mu_g  = gas_viscosity_air(T_film_K)
    k_g   = gas_thermal_conductivity_air(T_film_K)
    cp_g  = gas_cp_air(T_film_K)
    Dv    = gas_diffusivity_water_air(T_film_K)

    rho_l = liquid_density_water(Ts_K)
    cp_l  = liquid_cp_water(Ts_K)
    k_l   = liquid_thermal_conductivity_water(Ts_K)
    mu_l  = liquid_viscosity_water(Ts_K)
    alpha_l = k_l / (rho_l * cp_l)
    L_vap = latent_heat_water(Ts_K)

    # Solute / surfactant activity
    local x_w = 1.0
    if p.with_surfactants
        m_w = max(y[N+2], 1e-12)
        # Mole fraction Raoult
        n_w = m_w / M_WATER
        n_s = p.m_s_init / 0.05 # assume 50 g/mol typical surfactant
        x_w = n_w / (n_w + n_s)
    end

    # Sherwood & Nusselt numbers
    Sh = p.with_natural_convection ? sherwood_natural_convection(R, Ts_K, T_air_K, rho_g, mu_g, Dv) : 2.0
    Pr = calc_pr(mu_g, cp_g, k_g)
    Nu = p.with_natural_convection ? (2.0^3 + (0.54 * (max(GRAVITY * (1.0 / T_film_K) * abs(T_air_K - Ts_K) * ((2*R)^3) / ((mu_g/rho_g)^2) * Pr, 0.0)^0.25))^3)^(1.0/3.0) : 2.0

    # Evaporation rate
    P_sat_s = psat_water(Ts_K)
    Y_vs = surface_mass_fraction_from_psat(x_w * P_sat_s)
    P_sat_air = psat_water(T_air_K)
    Y_v_inf = surface_mass_fraction_from_psat(p.RH * P_sat_air)

    B_M = (Y_vs - Y_v_inf) / max(1.0 - Y_vs, 1e-12)
    bm_clamped = max(B_M, -0.99999)

    # Abramzon & Sirignano (1989) Eq. 11: m_dot = 2π * R * rho_g * Dv * Sh * ln(1 + B_M)
    m_dot = 2.0 * π * R * rho_g * Dv * Sh * log1p(bm_clamped)

    # Marangoni / surfactant factor
    if p.with_surfactants
        Gamma = y[N+3]
        Ma_norm = clamp(Gamma / p.Gamma_inf, 0.0, 1.0)
        m_dot *= (1.0 + 0.2 * (1.0 - Ma_norm))
    elseif p.with_marangoni
        m_dot *= marangoni_factor(R, Ts_K, T_air_K, mu_l, alpha_l)
    end

    # Interface retraction speed dR/dt
    dRdt = -m_dot / (4.0 * π * (R^2) * rho_l)
    dy[1] = dRdt

    # Heat conduction on immobilized Landau coordinate ξ ∈ [0, 1]
    dxi = 1.0 / (N - 1)
    dTdt = @view dy[2:N+1]

    # Center boundary condition (symmetry ∂T/∂ξ = 0, ∇²T = 3 ∂²T/∂ξ²)
    dTdt[1] = (alpha_l / (R^2)) * 6.0 * (T_nodes[2] - T_nodes[1]) / (dxi^2)

    # Interior nodes ξ_i = (i-1) * dxi
    inv_dxi2 = 1.0 / (dxi^2)
    inv_2dxi = 1.0 / (2.0 * dxi)
    inv_R2   = 1.0 / (R^2)
    dRdt_over_R = dRdt / R

    for i in 2:N-1
        xi = (i - 1) * dxi
        d2 = (T_nodes[i+1] - 2.0 * T_nodes[i] + T_nodes[i-1]) * inv_dxi2
        d1 = (T_nodes[i+1] - T_nodes[i-1]) * inv_2dxi
        lap = d2 + (2.0 / xi) * d1
        conv = dRdt_over_R * xi * d1
        dTdt[i] = alpha_l * inv_R2 * lap + conv
    end

    # Surface energy balance at ξ = 1
    # Convective heat transfer coefficient h_conv computed directly from Nusselt correlation
    h_conv = Nu * k_g / (2.0 * R)
    q_cond = k_l * (1.0 / R) * (T_nodes[N-1] - T_nodes[N]) / dxi
    q_conv = h_conv * (p.T_air - Ts_C)
    q_evap = m_dot * L_vap / (4.0 * π * (R^2))

    net_flux = q_cond + q_conv - q_evap
    dTdt[N] = net_flux / (rho_l * cp_l * (0.5 * R * dxi))

    # Surfactants / mass dynamics
    if p.with_surfactants
        dy[N+2] = -m_dot
        vol = (4.0 / 3.0) * π * (R^3)
        # Surfactant bulk molar concentration [mol/m^3] with M_surf = 0.05 kg/mol
        C_bulk = (p.m_s_init / 0.05) / vol
        Gamma = y[N+3]
        dG_dt = p.k_ads * C_bulk * (1.0 - Gamma / p.Gamma_inf) - p.k_des * Gamma
        dy[N+3] = dG_dt
    end

    return nothing
end

"""
    solve_radial_pde(cfg::RadialPDEConfig; solver=FBDF(), rtol=1e-6, atol=1e-8) -> RadialPDESolution

Integrate the 1D moving-boundary radial conduction PDE.
"""
function solve_radial_pde(cfg::RadialPDEConfig; solver=FBDF(), rtol=1e-6, atol=1e-8)
    N = cfg.n_nodes
    R0 = cfg.d0 / 2.0
    T0 = fill(cfg.T_drop_0, N)

    y0 = [R0; T0]
    if cfg.with_surfactants
        push!(y0, cfg.m_w_init)
        push!(y0, 0.0) # Gamma(0) = 0
    end

    cb = ContinuousCallback((u, t, integrator) -> u[1] - 1e-7, terminate!)

    prob = ODEProblem(radial_pde_rhs!, y0, (0.0, cfg.t_end), cfg)
    sol = solve(prob, solver; callback=cb, reltol=rtol, abstol=atol)

    n_steps = length(sol.t)
    R_vec = [sol.u[i][1] for i in 1:n_steps]
    T_mat = zeros(n_steps, N)
    for i in 1:n_steps
        T_mat[i, :] = sol.u[i][2:N+1]
    end

    m_w_vec = cfg.with_surfactants ? [sol.u[i][N+2] for i in 1:n_steps] : Float64[]
    Gamma_vec = cfg.with_surfactants ? [sol.u[i][N+3] for i in 1:n_steps] : Float64[]

    return RadialPDESolution(sol.t, R_vec, T_mat, m_w_vec, Gamma_vec)
end
