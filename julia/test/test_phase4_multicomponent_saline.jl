using Test
using Evaporacion
using LinearAlgebra

@testset "Phase 4: Multicomponent, Saline, Suspension & Crust Evaporation" begin

    # 1. Tonini & Cossali (2015) Stefan-Maxwell Exact Analytical Solution
    @testset "Stefan-Maxwell Analytical Matrix Formulation" begin
        # For single component (n=1), Tonini matches scalar Stefan blowing:
        # exp(N_hat) * (1 - y_s) = 1 - y_inf  =>  N_hat = ln((1 - y_inf) / (1 - y_s))
        y_s_val = 0.3
        y_inf_val = 0.05
        N_exact = log((1.0 - y_inf_val) / (1.0 - y_s_val))

        mu = [1.0]
        omega_species = reshape([1.0], 1, 1)
        omega_carrier = [1.0]
        W_s = [y_s_val]
        W_inf = [y_inf_val]

        res = stefan_maxwell_residual(N_exact, mu, W_s, W_inf, omega_species, omega_carrier)
        @test isapprox(res[1], 0.0; atol=1e-10)

        # Binary volatile mixture (n=2)
        mu2 = [0.6, 0.4]
        omega_sp2 = [1.0 1.2; 1.2 1.0]
        omega_c2 = [1.1, 0.9]
        A, B = stefan_maxwell_matrices(0.5, mu2, omega_sp2, omega_c2)
        @test size(A) == (2, 2)
        @test length(B) == 2
        @test A[1, 1] < 0.0
        @test A[2, 2] < 0.0
    end

    # 2. Ordoubadi (2019) Stefan-Fuchs Radial Convection Fluxes
    @testset "Stefan-Fuchs Multicomponent Fluxes" begin
        d = 30e-6 # 30 um
        rho_t = 1.18 # kg/m^3
        D_inert = 2.4e-5
        y_s_inert = 0.70
        y_inf_inert = 0.98

        # Two evaporating solvents (e.g. ethanol and water)
        y_s = [0.20, 0.10]
        y_inf = [0.01, 0.01]
        D_species = [1.2e-5, 2.5e-5]

        m_dot_t, m_dot_i = stefan_fuchs_multicomponent_fluxes(
            d, rho_t, D_inert, y_s_inert, y_inf_inert, y_s, y_inf, D_species
        )

        @test m_dot_t > 0.0
        @test length(m_dot_i) == 2
        @test all(m_dot_i .> 0.0)
        # Sum of individual fluxes is physically close to total (within ~25% for 2:1 diffusivity disparity)
        @test isapprox(sum(m_dot_i), m_dot_t; rtol=0.30)
    end

    # 3. Vehring-Ferron-Lechuga (VFL) Solute Surface Enrichment
    @testset "VFL Solute Surface Enrichment" begin
        # Zero evaporation rate -> no enrichment
        @test vfl_surface_enrichment(0.0) ≈ 1.0

        # Pe = 10 -> significant surface enrichment
        # E = 1 + 10/5 + 100/100 - 1000/4000 = 1 + 2 + 1 - 0.25 = 3.75
        E10 = vfl_surface_enrichment(10.0)
        @test isapprox(E10, 3.75; atol=1e-6)

        # Monotonic increase with Péclet number
        @test vfl_surface_enrichment(5.0) < vfl_surface_enrichment(10.0)
    end

    # 4. Fractal Nanoparticle Agglomerate Packing & Locking Diameter (Gerken 2014)
    @testset "Nanoparticle Packing & Crust Locking" begin
        # Gerken 2014 values: r_np = 25 nm, R_g = 125 nm, d_f = 1.7
        phi_max = fractal_agglomerate_max_packing(25e-9, 125e-9, 1.7)
        # Yields approx 0.0646 (6.46 vol%)
        @test isapprox(phi_max, 0.0646; atol=0.005)

        # Crust locking diameter
        d0 = 100e-6
        # Solid volume fraction 12.5% -> d_lock = d0 * (0.125)^(1/3) = 0.5 * d0
        d_lock = crust_locking_diameter(d0, 0.125; porosity=0.0)
        @test isapprox(d_lock, 50e-6; atol=1e-8)

        # Porous crust increases locking diameter
        d_lock_porous = crust_locking_diameter(d0, 0.125; porosity=0.5)
        @test d_lock_porous > d_lock
    end

    # 5. Multicomponent Droplet Simulation with Solute Enrichment & Crust
    @testset "Solve Multicomponent Droplet Simulation" begin
        water = MulticomponentSpecies(
            name = "water",
            molar_mass = 0.018015,
            density = 997.0,
            boiling_temp = 373.15,
            latent_heat = 2.26e6,
            diffusivity = 2.5e-5,
            is_volatile = true
        )
        solute = MulticomponentSpecies(
            name = "trehalose",
            molar_mass = 0.3423,
            density = 1580.0,
            is_volatile = false,
            liquid_diffusivity = 4e-10
        )

        cfg = MulticomponentDropletConfig(
            r0 = 20e-6, # 20 um radius
            T_d = 298.15,
            T_inf = 313.15,
            species = [water, solute],
            initial_mass_fractions = [0.90, 0.10]
        )

        sol = solve_multicomponent_droplet(cfg; t_span=(0.0, 0.1))
        @test sol.lifetime > 0.0
        @test sol.r[1] ≈ cfg.r0
        # Droplet radius decreases
        @test sol.r[end] < sol.r[1]
        # Solute mass fraction increases
        @test sol.mass_fractions[end, 2] > sol.mass_fractions[1, 2]
        # Solute surface enrichment is tracked and >= 1.0
        @test all(sol.surface_enrichments[:, 2] .>= 1.0)
        # Locking diameter is positive
        @test sol.locking_diameter > 0.0
    end

    # 6. Maruyama & Hasegawa (2020) Saline Vapor Pressure & Evaporation
    @testset "Saline Vapor Pressure & Evaporation Rate Constant" begin
        T_sur = 293.15
        T_inf = 293.15
        RH = 40.0

        # Pure water
        p_pure = saline_vapor_pressure(T_sur, 0.0)
        @test isapprox(p_pure, psat_water(T_sur); rtol=1e-5)

        # 10 wt% NaCl solution
        p_saline_10 = saline_vapor_pressure(T_sur, 0.10)
        @test p_saline_10 < p_pure

        # Evaporation constant beta decreases as salt concentrates
        beta_0 = saline_evaporation_constant(T_sur, T_inf, RH, 0.0)
        beta_10 = saline_evaporation_constant(T_sur, T_inf, RH, 0.10)
        beta_20 = saline_evaporation_constant(T_sur, T_inf, RH, 0.20)

        @test beta_0 > beta_10
        @test beta_10 > beta_20

        # Dry salt crystal diameter from mass conservation
        d0 = 1e-3 # 1 mm
        d_salt = saline_dry_crystal_diameter(d0, 0.10; rho_sol=1070.0, rho_salt=2160.0)
        @test d_salt > 0.0
        @test d_salt < d0
    end

    # 7. Transient Saline Droplet Simulation & Salt Precipitation
    @testset "Solve Saline Droplet Simulation" begin
        cfg_saline = SalineDropletConfig(
            d0 = 0.5e-3, # 0.5 mm
            T_sur = 295.15,
            T_inf = 298.15,
            RH = 30.0,
            Y_s0 = 0.10,
            solubility_limit = 0.2647
        )

        sol_saline = solve_saline_droplet(cfg_saline; t_span=(0.0, 150.0))
        @test sol_saline.precipitation_time > 0.0
        @test isapprox(sol_saline.d[1], cfg_saline.d0; atol=1e-8)
        @test sol_saline.d[end] < sol_saline.d[1]
        # Salt mass fraction increases to saturation limit
        @test sol_saline.Y_s[end] >= 0.25
        @test sol_saline.dry_salt_diameter > 0.0
        @test sol_saline.dry_salt_diameter < cfg_saline.d0
    end
end
