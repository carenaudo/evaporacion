"""
    Saline droplet evaporation and salt precipitation solvers.
    Based on Maruyama & Hasegawa (2020).
"""

using ..Common: R_UNIVERSAL, psat_water
using OrdinaryDiffEq

"""
    saline_vapor_pressure(T_s, Y_s; i=2.0, M_w=0.018015, M_s=0.05844)

Calculate the equilibrium vapor pressure [Pa] above an aqueous salt solution
depressed by Raoult's law with van 't Hoff electrolyte dissociation (Maruyama & Hasegawa 2020):
```math
Z_s = \\frac{Y_s M_w}{Y_s M_w + (1 - Y_s) M_s}
```
```math
P_{v,s} = P_{\\text{sat}}(T_s) (1 - i Z_s)
```
"""
function saline_vapor_pressure(T_s::Real, Y_s::Real;
                               i::Real=2.0,
                               M_w::Real=0.018015,
                               M_s::Real=0.05844)
    ys = clamp(float(Y_s), 0.0, 0.99)
    # Solute mole fraction
    Zs = (ys * float(M_w)) / (ys * float(M_w) + (1.0 - ys) * float(M_s))
    activity = max(0.0, 1.0 - float(i) * Zs)
    return psat_water(T_s) * activity
end

"""
    saline_evaporation_constant(T_sur, T_inf, RH, Y_s; D=2.5e-5, rho_l=997.0, i=2.0, M_w=0.018015, M_s=0.05844)

Compute the instantaneous evaporation constant ``\\beta(t) = -\\mathrm{d}(d^2)/\\mathrm{d}t`` [m^2/s]
from Maruyama & Hasegawa (2020) Eq. (7):
```math
\\beta(t) = \\frac{8 \\mathcal{D} M_w}{\\rho_l R_u} \\left[ \\frac{P_{\\text{sat}}(T_{\\text{sur}})(1 - i Z_s)}{T_{\\text{sur}}} - \\frac{P_{\\text{sat}}(T_{\\infty})}{T_{\\infty}} \\left(\\frac{RH}{100}\\right) \\right]
```
"""
function saline_evaporation_constant(T_sur::Real, T_inf::Real, RH::Real, Y_s::Real;
                                    D::Real=2.5e-5,
                                    rho_l::Real=997.0,
                                    i::Real=2.0,
                                    M_w::Real=0.018015,
                                    M_s::Real=0.05844)
    P_sur_eff = saline_vapor_pressure(T_sur, Y_s; i=i, M_w=M_w, M_s=M_s)
    P_inf_sat = psat_water(T_inf)
    P_inf_eff = P_inf_sat * (float(RH) / 100.0)

    driving_force = (P_sur_eff / float(T_sur)) - (P_inf_eff / float(T_inf))
    if driving_force <= 0.0
        return 0.0
    end

    return (8.0 * float(D) * float(M_w) / (float(rho_l) * R_UNIVERSAL)) * driving_force
end

"""
    saline_dry_crystal_diameter(d0, Y_s0; rho_sol=1070.0, rho_salt=2160.0)

Calculate the equivalent diameter of the precipitated solid salt crystal:
```math
d_{\\text{salt}} = d_0 \\left( \\frac{Y_{s,0} \\rho_{\\text{sol}}}{\\rho_{\\text{salt}}} \\right)^{1/3}
```
"""
function saline_dry_crystal_diameter(d0::Real, Y_s0::Real;
                                     rho_sol::Real=1070.0,
                                     rho_salt::Real=2160.0)
    vol_ratio = (float(Y_s0) * float(rho_sol)) / float(rho_salt)
    return float(d0) * cbrt(max(0.0, vol_ratio))
end

"""
    solve_saline_droplet(cfg::SalineDropletConfig{T}; t_span=nothing)

Solve transient saline droplet evaporation, solute concentration, and salt precipitation.
Integrates until the saturation solubility limit is exceeded, locking the shell diameter.
"""
function solve_saline_droplet(cfg::SalineDropletConfig{T};
                              t_span=nothing) where {T<:Real}
    # Initial solution density from ideal volume addition
    rho_sol0 = 1.0 / ((1.0 - cfg.Y_s0) / cfg.rho_water + cfg.Y_s0 / cfg.rho_salt)
    V0 = (pi / 6.0) * float(cfg.d0)^3
    m_tot0 = rho_sol0 * V0
    m_salt0 = cfg.Y_s0 * m_tot0
    m_water0 = m_tot0 - m_salt0

    d_salt = saline_dry_crystal_diameter(cfg.d0, cfg.Y_s0; rho_sol=rho_sol0, rho_salt=cfg.rho_salt)

    beta_init = saline_evaporation_constant(cfg.T_sur, cfg.T_inf, cfg.RH, cfg.Y_s0;
                                           D=cfg.diffusivity, rho_l=cfg.rho_water,
                                           i=cfg.i_van_t_hoff, M_w=cfg.water_molar_mass, M_s=cfg.salt_molar_mass)
    t_est = beta_init > 0.0 ? (float(cfg.d0)^2 / beta_init) : 1000.0
    t_end = t_span === nothing ? (1.5 * t_est) : t_span[2]

    # State: u = [m_water]
    function dm_dt(m_w, p, t)
        if m_w <= 1e-15
            return 0.0
        end
        # Instantaneous salt mass fraction
        ys = m_salt0 / (m_salt0 + m_w)
        # If precipitation reached, evaporation through crystal shell is locked or drops
        if ys >= cfg.solubility_limit
            return 0.0
        end

        V_curr = (m_w / cfg.rho_water) + (m_salt0 / cfg.rho_salt)
        d_curr = cbrt((6.0 * V_curr) / pi)

        beta = saline_evaporation_constant(cfg.T_sur, cfg.T_inf, cfg.RH, ys;
                                          D=cfg.diffusivity, rho_l=cfg.rho_water,
                                          i=cfg.i_van_t_hoff, M_w=cfg.water_molar_mass, M_s=cfg.salt_molar_mass)
        # d(m_w)/dt = rho_water * (pi/4) * d * (-beta)
        return -cfg.rho_water * (pi / 4.0) * d_curr * beta
    end

    condition(u, t, integrator) = (m_salt0 / (m_salt0 + u)) - cfg.solubility_limit
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    prob = ODEProblem(dm_dt, m_water0, (0.0, float(t_end)))
    sol = solve(prob, Tsit5(); callback=cb, reltol=1e-5, abstol=1e-9)

    t_vec = sol.t
    nt = length(t_vec)
    d_vec = zeros(Float64, nt)
    ys_vec = zeros(Float64, nt)
    beta_vec = zeros(Float64, nt)

    for step in 1:nt
        mw = max(1e-16, sol.u[step])
        ys_vec[step] = m_salt0 / (m_salt0 + mw)
        V_step = (mw / cfg.rho_water) + (m_salt0 / cfg.rho_salt)
        d_vec[step] = cbrt((6.0 * V_step) / pi)
        beta_vec[step] = saline_evaporation_constant(cfg.T_sur, cfg.T_inf, cfg.RH, ys_vec[step];
                                                    D=cfg.diffusivity, rho_l=cfg.rho_water,
                                                    i=cfg.i_van_t_hoff, M_w=cfg.water_molar_mass, M_s=cfg.salt_molar_mass)
    end

    precip_time = t_vec[end]

    return SalineSolution{T}(
        t_vec,
        d_vec,
        ys_vec,
        beta_vec,
        precip_time,
        d_salt
    )
end
