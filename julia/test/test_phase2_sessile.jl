using Test
using Evaporacion
using Evaporacion.Common
using Evaporacion.Sessile

@testset "Phase 2: Sessile Droplet Evaporation, Lifetimes, Spreading & Shielding" begin

    @testset "Picknett & Bexon (1977) Capacitance & Fluxes" begin
        # Small angle limit θ -> 0: f(θ) -> 2/π ≈ 0.6366
        f_0 = picknett_bexon_capacitance(0.0)
        @test isapprox(f_0, 2.0 / π, atol=1e-4)

        # Right angle θ = π/2: f(π/2) ≈ 1.0
        f_pi2 = picknett_bexon_capacitance(0.5 * π)
        @test isapprox(f_pi2, 1.0, atol=0.01)

        # Isolated sphere limit θ = π: f(π) = 4*ln(2)/π ≈ 0.8825
        f_pi = picknett_bexon_capacitance(π)
        @test isapprox(f_pi, 4.0 * log(2.0) / π, atol=0.01)

        # Bourges-Monnier & Shanahan (1995) flux
        R0 = 1e-3
        D = 2.6e-5
        delta_c = 0.02
        f_bms = bourges_monnier_flux(R0, 0.25 * π, D, delta_c)
        f_pb  = picknett_bexon_flux(R0, 0.25 * π, D, delta_c)
        @test f_bms > 0.0
        @test isapprox(f_bms, f_pb, rtol=0.15) # BMS is within 10-15% of PB
    end

    @testset "Thin-Droplet Analytical Lifetimes (Wilson 2023 Table 2)" begin
        R0 = 1e-3
        theta0 = 0.05 # small contact angle
        D = 2.6e-5
        delta_c = 0.02
        rho_l = 997.0

        t_cr_thin = lifetime_cr_thin(R0, theta0, D, delta_c, rho_l)
        t_ca_thin = lifetime_ca_thin(R0, theta0, D, delta_c, rho_l)
        @test t_cr_thin > 0.0
        @test isapprox(t_ca_thin / t_cr_thin, 1.5, atol=1e-8)

        # Compare with exact numerical integral lifetime_cr
        cfg_thin = SessileDropletConfig(D=D, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=theta0)
        t_cr_exact = lifetime_cr(cfg_thin)
        @test isapprox(t_cr_thin, t_cr_exact, rtol=0.01)

        # Stick-Slide lifetime ratio: (2θ₀ + θ*) / (2θ₀) = 1 + θ*/(2θ₀)
        theta_star = 0.02
        t_ss_thin = lifetime_ss_thin(R0, theta0, theta_star, D, delta_c, rho_l)
        expected_ratio = (2.0 * theta0 + theta_star) / (2.0 * theta0)
        @test isapprox(t_ss_thin / t_cr_thin, expected_ratio, atol=1e-8)

        # Stick-Jump lifetime
        t_sj_1 = lifetime_sj_thin(R0, theta0, theta_star, 1, D, delta_c, rho_l)
        @test isapprox(t_sj_1, t_cr_thin, atol=1e-8)
        t_sj_5 = lifetime_sj_thin(R0, theta0, theta_star, 5, D, delta_c, rho_l)
        @test t_sj_5 > t_cr_thin
    end

    @testset "Stauber Critical Angle & Hydrophobic Sphere" begin
        crit = stauber_critical_angle()
        @test isapprox(crit * 180.0 / π, 148.2, atol=0.1)

        V0 = 1e-9 # 1 mm^3
        t_pi = lifetime_hydrophobic_sphere(V0, 2.6e-5, 0.02, 997.0)
        @test t_pi > 0.0
    end

    @testset "Cooperative Shielding (Wray 2020 & Pair Lifetimes)" begin
        R0 = 1e-3
        theta0 = 0.05
        D = 2.6e-5
        delta_c = 0.02
        rho_l = 997.0

        # Exact pair lifetime at touching distance B = 2*R0:
        # arcsin(R0 / (2R0)) = arcsin(1/2) = π/6
        # t_pair = [1 + (2/π)*(π/6)] * t_CR,thin = (1 + 1/3) * t_CR,thin = 4/3 * t_CR,thin
        t_pair_touch = lifetime_pair_cr(2.0 * R0, R0, theta0, D, delta_c, rho_l)
        t_cr = lifetime_cr_thin(R0, theta0, D, delta_c, rho_l)
        @test isapprox(t_pair_touch / t_cr, 4.0 / 3.0, atol=1e-7)

        # Cooperative shielding linear system for 2 identical droplets at B = 2*R0
        radii = [R0, R0]
        positions = [0.0 0.0; 2.0 * R0 0.0]
        thetas = [theta0, theta0]
        fluxes = cooperative_shielding_fluxes(radii, positions, D, delta_c, thetas)

        F_iso = popov_global_flux(R0, theta0, D, delta_c)
        # Each droplet flux must be exactly 0.75 * F_iso because M = [1 1/3; 1/3 1] => F = [3/4; 3/4] * F_iso
        @test isapprox(fluxes[1] / F_iso, 0.75, atol=1e-7)
        @test isapprox(fluxes[2] / F_iso, 0.75, atol=1e-7)
    end

    @testset "Shallow Wells & Conduction-Limited Regime" begin
        # Well lifetime
        R0 = 1e-3
        theta0 = 0.1
        t_well = lifetime_well_cr(0.5e-4, 2.0, R0, theta0, 2.6e-5, 0.02, 997.0)
        t_flat = lifetime_cr_thin(R0, theta0, 2.6e-5, 0.02, 997.0)
        @test t_well > t_flat

        # Conduction-limited lifetime
        t_cond = lifetime_conduction_limited(R0, theta0, 997.0, 2.4e6, 1e-3, 1.4, 0.02, 0.001)
        @test t_cond > 0.0
    end

    @testset "Saxton (2016) Spreading-Evaporation & d^(13/7) Law" begin
        # Effective contact angle Φ₀(E)
        @test isapprox(saxton_effective_contact_angle(0.0), 1.0, atol=1e-6)
        phi_1 = saxton_effective_contact_angle(1.0, 1)
        @test phi_1 > 1.0

        # Near extinction scaling S₀(T)
        Tc = 10.0
        S_9 = saxton_extinction_scaling(9.0, Tc, 1.0)
        S_99 = saxton_extinction_scaling(9.9, Tc, 1.0)
        @test S_9 > S_99
        @test isapprox(saxton_extinction_scaling(10.0, Tc, 1.0), 0.0, atol=1e-10)

        # Spreading ODE solver
        sol_sax = solve_saxton_spreading_evaporation(1.0, 0.1, (0.0, 2.0); s0_init=1.0)
        @test Symbol(sol_sax.retcode) == :Success
        @test sol_sax.u[end][1] > 0.0
    end

    @testset "Unpinned CCA Flow & Marangoni Inversion Diagnostics" begin
        # CCA depth-averaged velocity
        R0 = 1e-3
        theta0 = 0.5
        u_cca = depth_averaged_radial_velocity_cca(0.5 * R0, R0, theta0, 2.6e-5, 0.02, 997.0)
        @test u_cca > 0.0

        # Rocha et al. (2024) Marangoni inversion
        # Clean water: Gamma0 = 0 -> no flow inversion
        mar_clean = marangoni_thermal_solutal_analysis(R0, 0.8, 0.0, 1e-3, 1.4e-7, -1.6e-4, -100.0)
        @test !mar_clean.flow_inverted

        # Trace contamination: solutal Marangoni overwhelms thermal Marangoni
        mar_contam = marangoni_thermal_solutal_analysis(R0, 0.8, 1e-5, 1e-3, 1.4e-7, -1.6e-4, -5e4)
        @test mar_contam.flow_inverted
    end

end
