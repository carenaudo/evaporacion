"""
    Transient multicomponent droplet evaporation solver.
"""

using ..Common: R_UNIVERSAL, psat_water, gas_density_air
using OrdinaryDiffEq

"""
    solve_multicomponent_droplet(cfg::MulticomponentDropletConfig{T}; t_span=nothing, dt_save=1e-4)

Solve transient multicomponent droplet evaporation and solute surface enrichment.
"""
function solve_multicomponent_droplet(cfg::MulticomponentDropletConfig{T};
                                      t_span=nothing,
                                      dt_save::Real=1e-4) where {T<:Real}
    n_sp = length(cfg.species)
    # Initial droplet volume and mass
    V0 = (4.0 / 3.0) * pi * float(cfg.r0)^3
    # Initial component masses from ideal volume addition
    rho_l_mix = 1.0 / sum(cfg.initial_mass_fractions[i] / cfg.species[i].density for i in 1:n_sp)
    m_tot0 = rho_l_mix * V0
    m0_vec = [float(cfg.initial_mass_fractions[i]) * m_tot0 for i in 1:n_sp]

    # Non-volatile volume fraction
    V_solute0 = sum(cfg.species[i].is_volatile ? 0.0 : m0_vec[i] / cfg.species[i].density for i in 1:n_sp)
    phi_solute0 = V0 > 0.0 ? V_solute0 / V0 : 0.0
    d_lock = crust_locking_diameter(2.0 * cfg.r0, phi_solute0)

    # Carrier gas density
    rho_g = (float(cfg.P_tot) * cfg.carrier_molar_mass) / (R_UNIVERSAL * float(cfg.T_inf))

    t_end = t_span === nothing ? 1.0 : t_span[2]

    # ODE state: u = [m_1, m_2, ..., m_n]
    function dm_dt!(du, u, p, t)
        m_tot = sum(u)
        if m_tot <= 1e-18
            du .= 0.0
            return
        end

        # Liquid volume and instantaneous radius
        V_liq = sum(u[i] / cfg.species[i].density for i in 1:n_sp)
        r_inst = cbrt((3.0 * V_liq) / (4.0 * pi))
        d_inst = 2.0 * r_inst

        # Check if crust has locked the droplet diameter
        if d_inst <= d_lock && phi_solute0 > 0.0
            # Volatile evaporation through crust proceeds at slower rate or stops
            du .= 0.0
            return
        end

        # Liquid mass fractions and mole fractions
        # Liquid mole fractions via Raoult's law for ideal liquid mixture (Tonini & Cossali 2015)
        moles_liq = [max(0.0, u[i]) / cfg.species[i].molar_mass for i in 1:n_sp]
        tot_moles = max(sum(moles_liq), 1e-15)
        X_liq = [moles_liq[i] / tot_moles for i in 1:n_sp]

        # Surface vapor partial pressures via Raoult's law for volatile components
        # p_i = X_liq[i] * p_sat_i
        p_vap = zeros(Float64, n_sp)
        for i in 1:n_sp
            if cfg.species[i].is_volatile
                # Approximate psat via Antoine or water psat / boiling temp
                p_sat_i = if cfg.species[i].name == "water"
                    psat_water(cfg.T_d)
                else
                    # Clausius-Clapeyron from boiling point
                    P_ATM_CONST = 101325.0
                    R_spec = R_UNIVERSAL / cfg.species[i].molar_mass
                    P_ATM_CONST * exp((cfg.species[i].latent_heat / R_spec) * (1.0 / cfg.species[i].boiling_temp - 1.0 / cfg.T_d))
                end
                p_vap[i] = X_liq[i] * p_sat_i
            end
        end

        # Surface vapor mass fractions in gas phase
        P_tot = float(cfg.P_tot)
        p_inert = max(1.0, P_tot - sum(p_vap))
        m_inert_gas = p_inert * cfg.carrier_molar_mass
        m_vap_gas = sum(p_vap[i] * cfg.species[i].molar_mass for i in 1:n_sp)
        m_gas_tot = m_inert_gas + m_vap_gas

        y_s = [(p_vap[i] * cfg.species[i].molar_mass) / m_gas_tot for i in 1:n_sp]
        y_s_inert = m_inert_gas / m_gas_tot
        y_inf = if isempty(cfg.ambient_mole_fractions)
            zeros(Float64, n_sp)
        else
            x_inf = [float(cfg.ambient_mole_fractions[i]) for i in 1:n_sp]
            x_carrier = max(0.0, 1.0 - sum(x_inf))
            M_mix_inf = sum(x_inf[i] * cfg.species[i].molar_mass for i in 1:n_sp) + x_carrier * cfg.carrier_molar_mass
            [(x_inf[i] * cfg.species[i].molar_mass) / max(M_mix_inf, 1e-12) for i in 1:n_sp]
        end
        y_inf_inert = max(0.01, 1.0 - sum(y_inf))

        D_species = [cfg.species[i].diffusivity for i in 1:n_sp]
        _, m_dot_i = stefan_fuchs_multicomponent_fluxes(
            d_inst, rho_g, cfg.carrier_diffusivity,
            y_s_inert, y_inf_inert, y_s, y_inf, D_species
        )

        for i in 1:n_sp
            if cfg.species[i].is_volatile && u[i] > 1e-15
                du[i] = -m_dot_i[i]
            else
                du[i] = 0.0
            end
        end
    end

    condition(u, t, integrator) = sum(u) - 1e-16
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    prob = ODEProblem(dm_dt!, m0_vec, (0.0, float(t_end)))
    sol = solve(prob, Tsit5(); callback=cb, reltol=1e-5, abstol=1e-9)

    t_vec = sol.t
    nt = length(t_vec)
    masses = zeros(Float64, nt, n_sp)
    mass_fracs = zeros(Float64, nt, n_sp)
    r_vec = zeros(Float64, nt)
    enrichments = zeros(Float64, nt, n_sp)

    for step in 1:nt
        u_step = sol.u[step]
        m_tot_step = max(1e-18, sum(u_step))
        V_step = sum(u_step[i] / cfg.species[i].density for i in 1:n_sp)
        r_vec[step] = cbrt((3.0 * V_step) / (4.0 * pi))
        for i in 1:n_sp
            masses[step, i] = max(0.0, u_step[i])
            mass_fracs[step, i] = masses[step, i] / m_tot_step
        end

        # Calculate kappa and VFL surface enrichment for solutes
        kappa_step = if step > 1
            dt = t_vec[step] - t_vec[step-1]
            dt > 0.0 ? -( (2.0 * r_vec[step])^2 - (2.0 * r_vec[step-1])^2 ) / dt : 0.0
        else
            0.0
        end

        for i in 1:n_sp
            if !cfg.species[i].is_volatile
                pe = solute_peclet_number(kappa_step, cfg.species[i].liquid_diffusivity)
                enrichments[step, i] = vfl_surface_enrichment(pe)
            else
                enrichments[step, i] = 1.0
            end
        end
    end

    lifetime = t_vec[end]

    return MulticomponentSolution{T}(
        t_vec,
        r_vec,
        masses,
        mass_fracs,
        enrichments,
        lifetime,
        d_lock
    )
end
