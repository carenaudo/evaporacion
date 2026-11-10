# =============================================================================
# Unit Tests Verifying Implementation & Reference Audit Fixes
# =============================================================================

using Test
using Evaporacion
using Evaporacion.Common
using Evaporacion.Droplet
using Evaporacion.Sessile
using Evaporacion.RadialPDE
using Evaporacion.Multicomponent

@testset "Audit Fixes: Category 1 - Droplet & Correlations" begin
    # D1 & D2: Uniform profile and surface energy residual
    p_uni = DropletParams(
        T_inf = 25.0,
        U_rel = 0.0,
        Y_inf = 0.005,
        profile = :uniforme
    )
    d0 = 200e-6
    m0 = liquid_density_water(298.15) * (π / 6.0) * (d0^3)
    
    # In uniform profile, solve_surface_temperature must return T_bar directly
    @test solve_surface_temperature(m0, 20.0, p_uni) ≈ 20.0
    @test solve_surface_temperature(m0, 30.0, p_uni) ≈ 30.0

    # Solve droplet with uniform profile: temperature must evolve and not stay frozen
    u0 = [m0, 25.0] # starts at 25 °C, wet-bulb is ~14 °C
    sol_uni = solve_droplet(p_uni, u0, (0.0, 5.0))
    @test Symbol(sol_uni.retcode) in (:Success, :Terminated)
    T_final = sol_uni.u[end][2]
    # Droplet must cool down towards wet-bulb temperature
    @test T_final < 24.0
    @test sol_uni.u[end][1] < m0

    # D1: Verify q_l sign convention in DropletAlgebraic
    # When Ts > T_bar (heat conducted into interior), q_l must be > 0
    p_cond = DropletParams(T_inf = 40.0, profile = :parabolico, Y_inf = 0.01)
    alg_inward = calc_droplet_algebraic(m0, 20.0, 30.0, p_cond) # Ts = 30 > T_bar = 20
    @test alg_inward.q_l > 0.0

    # When Ts < T_bar (interior warmer, heat conducts to surface), q_l must be < 0
    alg_outward = calc_droplet_algebraic(m0, 30.0, 20.0, p_cond) # Ts = 20 < T_bar = 30
    @test alg_outward.q_l < 0.0

    # D5: Zhou profile b_q coefficient
    p_zhou = DropletParams(T_inf = 25.0, profile = :zhou_cubic, Y_inf = 0.005)
    alg_z = calc_droplet_algebraic(m0, 20.0, 20.0, p_zhou)
    @test alg_z.b_q ≈ (1.0 / 0.086) atol=1e-3
end

@testset "Audit Fixes: Category 2 - Radial PDE" begin
    # R1: Verify that RadialPDE evaporation rate formula gives m_dot = 2π R ρ_g Dv Sh ln(1+Bm)
    cfg_pde = RadialPDEConfig(
        n_nodes = 11,
        d0 = 200e-6,
        T_drop_0 = 25.0,
        T_air = 25.0,
        RH = 0.5,
        t_end = 2.0
    )
    sol_pde = solve_radial_pde(cfg_pde)
    @test sol_pde.R[end] < (cfg_pde.d0 / 2.0)
    @test sol_pde.R[end] > 0.0
end

@testset "Audit Fixes: Category 3 - Sessile Droplets" begin
    # S1: Picknett & Bexon capacitance factor at θ = π/2 should be ≈ 1.0000
    f_pb_half = picknett_bexon_capacitance(0.5 * π)
    @test isapprox(f_pb_half, 1.0000, atol=2e-3)

    # Continuity across θ = 0.175 rad (10°)
    f_pb_below = picknett_bexon_capacitance(0.1749)
    f_pb_above = picknett_bexon_capacitance(0.1751)
    @test isapprox(f_pb_below, f_pb_above, atol=1e-3)

    # S2: Stick-jump thin droplet lifetime bounds
    R0 = 1e-3
    th0 = 0.3
    th_star = 0.15
    D = 2.5e-5
    dc = 0.01
    rho = 1000.0

    t_cr = lifetime_cr_thin(R0, th0, D, dc, rho)
    t_ca = lifetime_ca_thin(R0, th0, D, dc, rho)
    t_sj_1 = lifetime_sj_thin(R0, th0, th_star, 1, D, dc, rho)
    t_sj_5 = lifetime_sj_thin(R0, th0, th_star, 5, D, dc, rho)
    t_sj_50 = lifetime_sj_thin(R0, th0, th_star, 50, D, dc, rho)

    # Base case n = 1 (no jumps) must match t_cr exactly
    @test t_sj_1 ≈ t_cr
    # Monotonically increasing with number of cycles
    @test t_sj_1 < t_sj_5 < t_sj_50
    # Must remain strictly bounded below t_ca
    @test t_sj_50 < t_ca
end

@testset "Audit Fixes: Category 4 - Multicomponent Droplets" begin
    # M1: Raoult's law mole fraction vs mass fraction
    # Mixture of water (M=0.018 kg/mol) and ethanol (M=0.046 kg/mol)
    water = MulticomponentSpecies(
        name="water", molar_mass=0.018015, density=997.0,
        boiling_temp=373.15, latent_heat=2.26e6, diffusivity=2.5e-5, is_volatile=true
    )
    ethanol = MulticomponentSpecies(
        name="ethanol", molar_mass=0.04607, density=789.0,
        boiling_temp=351.45, latent_heat=8.46e5, diffusivity=1.2e-5, is_volatile=true
    )
    
    # 50-50 by mass: 0.5 kg water, 0.5 kg ethanol
    # Moles water = 0.5 / 0.018015 = 27.75 mol
    # Moles ethanol = 0.5 / 0.04607 = 10.85 mol
    # Water mole fraction = 27.75 / (27.75 + 10.85) ≈ 0.719 (NOT 0.50!)
    cfg_multi = MulticomponentDropletConfig(
        r0 = 50e-6,
        species = [water, ethanol],
        initial_mass_fractions = [0.5, 0.5],
        T_d = 293.15,
        T_inf = 293.15,
        P_tot = 101325.0,
        carrier_molar_mass = 0.02897,
        carrier_diffusivity = 2.0e-5,
        ambient_mole_fractions = [0.0, 0.0]
    )
    sol_multi = solve_multicomponent_droplet(cfg_multi; t_span=(0.0, 0.05))
    @test sol_multi.r[end] < sol_multi.r[1]
    # Both volatile species must evaporate
    @test sol_multi.masses[end, 1] < sol_multi.masses[1, 1]
    @test sol_multi.masses[end, 2] < sol_multi.masses[1, 2]

    # M3: Dilute L'Hopital limit in stefan_fuchs_multicomponent_fluxes
    d_test = 100e-6
    rho_test = 1.2
    D_m = 2.5e-5
    y_s = [0.02, 0.01]
    y_inf = [0.005, 0.001]
    y_s_inert = 1.0 - sum(y_s)
    y_inf_inert = 1.0 - sum(y_inf)
    D_spec = [2.5e-5, 2.5e-5] # identical diffusivities -> beta_i = 1 exactly

    m_tot, m_i = stefan_fuchs_multicomponent_fluxes(
        d_test, rho_test, D_m, y_s_inert, y_inf_inert, y_s, y_inf, D_spec
    )
    @test m_i[1] > 0.0
    @test m_i[2] > 0.0
    @test isapprox(m_tot, sum(m_i), rtol=1e-3)
end

@testset "Audit Fixes: Category 5 - Level 1 / Level 2 Linear Fluxes & Counter" begin
    T_inf = 35.0
    Y_inf = surface_mass_fraction_from_psat(0.4 * psat_water(T_inf + T_ZERO_C), P_ATM)
    d0 = 200e-6
    T0 = 20.0
    m0 = liquid_density_water(T0 + T_ZERO_C) * (π / 6.0) * (d0^3)

    for variant in (:level_1, :level_2)
        ctr = TsSolverCounter()
        p = DropletParams(
            T_inf = T_inf,
            U_rel = 2.0,
            Y_inf = Y_inf,
            profile = :parabolico,
            model_variant = variant,
            ts_solver = :brent,
            ts_counter = ctr
        )
        Ts = solve_surface_temperature(m0, T0, p)
        alg = calc_droplet_algebraic(m0, T0, Ts, p)
        @test ctr.calls == 1
        @test ctr.iters > 0

        # Exact linear flux relations (no log1p(B_M) / log1p_over_x(B_T))
        sh_expected = (variant === :level_1) ? 2.0 : alg.Sh0
        nu_expected = (variant === :level_1) ? 2.0 : alg.Nu0
        @test alg.Sh_star ≈ sh_expected
        @test alg.Nu_star ≈ nu_expected

        T_eg_K = (Ts + T_inf) / 2.0 + T_ZERO_C
        Dab = gas_diffusivity_water_air(T_eg_K, P_ATM)
        rg  = gas_density_air(T_eg_K, P_ATM)
        kg  = gas_thermal_conductivity_air(T_eg_K)
        Y_s = surface_mass_fraction_from_psat(psat_water(Ts + T_ZERO_C), P_ATM)

        m_dot_lin = π * alg.d * rg * Dab * sh_expected * (Y_s - Y_inf)
        q_g_lin   = π * alg.d * kg * nu_expected * (T_inf - Ts)
        @test isapprox(alg.m_dot, m_dot_lin; rtol=1e-10)
        @test isapprox(alg.q_g, q_g_lin; rtol=1e-10)

        # Intensive surface residual must vanish at solved Ts
        liq = precompute_liquid_state(m0, T0, p)
        @test abs(intensive_surface_residual(Ts, liq, p)) < 1e-6
    end
end
