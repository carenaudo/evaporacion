using Test
using Evaporacion
using Evaporacion.Common
using Evaporacion.Droplet

@testset "Phase 1: Spherical Evaporation & Thermal Closures" begin

    @testset "Thermal Expansion & 1/3 Reference State" begin
        # Water density maximum near 4 °C (277.15 K) -> β ≈ 0
        beta_4 = liquid_thermal_expansion_water(277.15)
        @test abs(beta_4) < 5e-5

        # Room temperature (20 °C = 293.15 K) -> β ≈ 2.07e-4 K^-1
        beta_20 = liquid_thermal_expansion_water(293.15)
        @test isapprox(beta_20, 2.07e-4, rtol=0.1)

        # 1/3 Reference state
        ref = reference_state_one_third(290.0, 320.0, 0.01, 0.002)
        @test isapprox(ref.T_ref, 300.0, atol=1e-6)
        @test isapprox(ref.Y_v_ref, 0.01 + (1.0 / 3.0) * (0.002 - 0.01), atol=1e-8)

        # Invertibility of mole <-> mass fractions
        for y in [0.001, 0.01, 0.05, 0.2]
            x = mole_fraction_from_mass_fraction(y)
            y_back = mass_fraction_from_mole_fraction(x)
            @test isapprox(y, y_back, atol=1e-7)
        end
    end

    @testset "Wilke Mixing & Gas Mixture Properties" begin
        # Pure dry air limit (Y_v = 0)
        mix_dry = gas_mixture_properties(293.15, 0.0)
        @test isapprox(mix_dry.mu_g, gas_viscosity_air(293.15), rtol=1e-5)
        @test isapprox(mix_dry.k_g, gas_thermal_conductivity_air(293.15), rtol=1e-5)
        @test isapprox(mix_dry.rho_g, gas_density_air(293.15), rtol=1e-5)
        @test isapprox(mix_dry.cp_g, gas_cp_air(293.15), rtol=1e-5)

        # Humid mixture properties at 40 °C (313.15 K)
        mix_humid = gas_mixture_properties(313.15, 0.03)
        @test mix_humid.rho_g > 0.9 && mix_humid.rho_g < 1.3
        @test mix_humid.mu_g > 1.5e-5 && mix_humid.mu_g < 2.5e-5
        @test mix_humid.k_g > 0.02 && mix_humid.k_g < 0.05
        @test mix_humid.D_ab > 2.0e-5
    end

    @testset "Spheroid Deformation & Holterman Correlations" begin
        # Spheroid deformation ratio at ε = 0 (sphere)
        @test isapprox(spheroid_deformation_evaporation_ratio(0.0), 1.0, atol=1e-5)

        # Oblate spheroid at ε = 0.5 (a/b = 1.5): Tonini & Cossali (2015) exact formula
        ratio_05 = spheroid_deformation_evaporation_ratio(0.5)
        @test isapprox(ratio_05, 1.161, atol=0.01)

        # Holterman evaporation rate & relaxation time
        d_rate = holterman_evaporation_rate(100e-6, 2.0, 1e-10, 0.552)
        @test d_rate < 0.0 # negative shrinkage rate

        tau_c = holterman_relaxation_time(100e-6, 1000.0, 4184.0, 0.026, 2.5e-5, 18.015, 2.45e6, 293.15, 66.0)
        @test tau_c > 0.0 && tau_c < 1.0
    end

    @testset "Ventilated Wet-Bulb & Revised D2 Law" begin
        # Ventilated wet-bulb temperature at 25 °C, Re=10
        T_ds = ventilated_wet_bulb_temperature(25.0, 0.005, 10.0)
        @test T_ds > 5.0 && T_ds < 25.0

        # Revised D2-law rate
        rev = calc_d2_revised_rate(25.0, 0.005, 10.0, 100e-6)
        @test rev.Kr > 0.0
        @test rev.tau > 0.0 && rev.tau < 100.0
        @test rev.Sh > 2.0
    end

    @testset "Zhou (2017) 3rd-Order Profile & Hill's Vortex Stagnation" begin
        rs = 1e-3
        T_bar = 20.0
        Ts = 30.0

        # Reconstruct profile
        T_0 = reconstruct_zhou_cubic_profile(0.0, rs, T_bar, Ts)
        T_s_rec = reconstruct_zhou_cubic_profile(rs, rs, T_bar, Ts)
        @test isapprox(T_s_rec, Ts, atol=1e-8)

        # Test zero derivative at Hill's vortex stagnation center: r = rs / sqrt(2)
        r_vortex = rs / sqrt(2.0)
        dr = 1e-7 * rs
        dT_dr = (reconstruct_zhou_cubic_profile(r_vortex + dr, rs, T_bar, Ts) -
                 reconstruct_zhou_cubic_profile(r_vortex - dr, rs, T_bar, Ts)) / (2.0 * dr)
        @test isapprox(dT_dr, 0.0, atol=1e-5)
    end

    @testset "Sazhin (2006, 2017) Analytical Transient Eigenvalue Series" begin
        rs = 1e-3
        T0 = 20.0
        T_eff = 80.0
        kappa_eff = 1.4e-7 # thermal diffusivity of water
        Bi = 2.0

        # Initial time t = 0 -> T ≈ T0
        T_init = sazhin_eigenvalue_temperature_profile(0.0, rs, 0.0, T0, T_eff, kappa_eff, Bi; n_modes=50)
        @test isapprox(T_init, T0, atol=1.0)

        # Long time t -> ∞ -> T -> T_eff
        T_long = sazhin_eigenvalue_temperature_profile(0.0, rs, 100.0, T0, T_eff, kappa_eff, Bi; n_modes=30)
        @test isapprox(T_long, T_eff, atol=1e-3)
    end

    @testset "ODE Solver with Phase 1 Features" begin
        # DropletParams with Zhou cubic profile, 1/3 reference state, Wilke mixture, thermal dilation
        params_p0 = DropletParams(
            T_inf = 40.0,
            U_rel = 1.5,
            Y_inf = 0.005,
            profile = :zhou_cubic,
            reference_state = :one_third,
            use_wilke_mixtures = true,
            thermal_dilation = true,
            variable_cp_factor = true,
            spheroid_epsilon = 0.1
        )

        d0 = 100e-6
        m0 = liquid_density_water(293.15) * (π / 6.0) * (d0^3)
        u0 = [m0, 20.0]
        tspan = (0.0, 0.5)

        sol = solve_droplet(params_p0, u0, tspan)
        ret_sym = Symbol(sol.retcode)
        @test ret_sym == :Success || ret_sym == :Terminated
        @test sol.u[end][1] < m0 # Mass evaporated
        @test sol.u[end][2] > 15.0 # Droplet temperature realistic
    end

end
