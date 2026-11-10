using Test
using Evaporacion

@testset "Phase 3: Kinetic & Knudsen-Layer Evaporation" begin
    # 1. Kelvin curvature shift
    @testset "Kelvin Curvature Shift" begin
        # Macroscopic droplet (r = 100 um) -> negligible Kelvin elevation
        k_macro = kelvin_vapor_pressure_ratio(100e-6, 293.15)
        @test isapprox(k_macro, 1.0; atol=1e-4)

        # Nanoscale droplet (r = 10 nm) -> significant vapor pressure elevation
        k_nano = kelvin_vapor_pressure_ratio(10e-9, 293.15)
        @test k_nano > 1.10
        @test k_nano < 1.15

        # Saturation vapor pressure elevates with curvature
        p_flat = psat_water(293.15)
        p_curv = psat_kelvin(10e-9, 293.15)
        @test p_curv ≈ p_flat * k_nano
    end

    # 2. Molecular thermal speed & Vlasov transition length
    @testset "Molecular Thermal Speed & Transition Length" begin
        v_th = mean_thermal_speed(293.15; molar_mass=0.018015)
        # For water at 20 C, mean molecular speed is ~586.8 m/s
        @test isapprox(v_th, 586.8; atol=2.0)

        D = 2.5e-5 # m^2/s
        alpha = 1.0
        lk = vlasov_kinetic_length(D, alpha, v_th)
        # l_k = 4D / (alpha * v) ~ 0.17 um
        @test isapprox(lk, 1.7e-7; rtol=0.1)

        # Reduced accommodation coefficient alpha increases kinetic length
        lk_small_alpha = vlasov_kinetic_length(D, 0.1, v_th)
        @test lk_small_alpha ≈ 10.0 * lk
    end

    # 3. Vlasov Evaporation Rate (Continuum vs Free Molecular Asymptotics)
    @testset "Vlasov Evaporation Rate Asymptotics" begin
        D = 2.5e-5
        delta_c = 0.5 # mol/m^3
        v_th = 587.0
        alpha = 1.0
        lk = vlasov_kinetic_length(D, alpha, v_th)

        # Continuum limit (r >> lk)
        r_large = 1e-3 # 1 mm >> 0.17 um
        I_large = vlasov_evaporation_rate(r_large, D, delta_c, alpha, v_th)
        I_maxwell = 4.0 * pi * r_large * D * delta_c
        @test isapprox(I_large, I_maxwell; rtol=1e-3)

        # Kinetic / Free-molecular limit (r << lk)
        r_small = 1e-9 # 1 nm << 0.17 um
        I_small = vlasov_evaporation_rate(r_small, D, delta_c, alpha, v_th)
        I_knudsen = pi * r_small^2 * alpha * v_th * delta_c
        @test isapprox(I_small, I_knudsen; rtol=1e-2)
    end

    # 4. Vlasov Analytical Trajectory & Lifetime Extension
    @testset "Vlasov Analytical Trajectory & Lifetime" begin
        r0 = 2e-6 # 2 um
        D = 2.4e-5
        delta_c = 0.4
        rho_l = 997.0
        M_v = 0.018015
        alpha = 1.0
        v_th = mean_thermal_speed(293.15)
        lk = vlasov_kinetic_length(D, alpha, v_th)

        # Lifetime
        tau_vlasov = vlasov_lifetime(r0, D, delta_c, rho_l, M_v, alpha, v_th)
        tau_maxwell = (rho_l * r0^2) / (2.0 * D * M_v * delta_c)

        # Kinetic resistance extends droplet lifetime
        @test tau_vlasov > tau_maxwell
        @test isapprox(tau_vlasov, tau_maxwell * (1.0 + 2.0 * lk / r0); rtol=1e-6)

        # At t = 0, R(0) = r0
        R_at_0 = vlasov_radius_trajectory(0.0, r0, D, delta_c, rho_l, M_v, alpha, v_th)
        @test isapprox(R_at_0, r0; rtol=1e-6)

        # At t = tau, R(tau) = 0
        R_at_tau = vlasov_radius_trajectory(tau_vlasov, r0, D, delta_c, rho_l, M_v, alpha, v_th)
        @test isapprox(R_at_tau, 0.0; atol=1e-12)
    end

    # 5. Zhao & Nadal (2023) Kinetic Péclet & Lifetime Extension
    @testset "Zhao & Nadal (2023) Kinetic Péclet & Extension" begin
        D = 2.4e-5
        r0 = 2e-6
        T_inf = 293.15
        Pe0 = kinetic_peclet_number(r0, T_inf, D)
        @test Pe0 > 0.0

        tau_d2 = 0.005 # s
        tau_ext = zhao_nadal_lifetime(tau_d2, Pe0)
        @test tau_ext > tau_d2
    end

    # 6. Akkuş (2020) Stefan Blowing Factor
    @testset "Akkuş (2020) Stefan Flow Factor" begin
        # Zero vapor mass fraction
        @test stefan_flow_enhancement(0.0) ≈ 1.0
        # Moderate vapor mass fraction
        @test isapprox(stefan_flow_enhancement(0.1), 1.1111; atol=1e-3)
        # High vapor mass fraction (e.g. heated droplet)
        @test stefan_flow_enhancement(0.5) ≈ 2.0
    end

    # 7. Transient Kinetic Droplet ODE Solution
    @testset "Solve Transient Kinetic Droplet ODE" begin
        cfg = KineticDropletConfig(
            r0 = 3e-6,
            T_s = 293.15,
            T_inf = 298.15,
            RH = 0.4,
            alpha = 1.0,
            use_kelvin = true,
            use_stefan = true
        )
        sol = solve_kinetic_droplet(cfg)

        @test sol.lifetime > 0.0
        @test sol.r[1] ≈ cfg.r0
        @test sol.r[end] < 1e-7
        # Droplet radius decreases monotonically
        @test all(diff(sol.r) .<= 0.0)
        # Non-equilibrium kinetic effects extend lifetime beyond continuum Maxwellian
        @test sol.lifetime >= sol.tau_continuum
    end
end
