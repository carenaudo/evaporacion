using Test
using Evaporacion

@testset "Phase 5: Evaporative Extinction & Atmospheric Coupling" begin

    # 1. Renaudo et al. (2022, 2023) Evaporation Parameter k
    @testset "Renaudo Evaporation Parameter" begin
        # 25 °C, 50% RH
        k_50 = renaudo_evaporation_parameter(25.0, 50.0)
        @test k_50 > 0.0

        # Drier air (20% RH) has larger wet-bulb depression, hence smaller k (faster evaporation)
        k_20 = renaudo_evaporation_parameter(25.0, 20.0)
        @test k_20 < k_50
    end

    # 2. Renaudo Cut-off Extinction Diameter d_min
    @testset "Renaudo Extinction Diameter d_min" begin
        H = 0.5 # 50 cm release height
        T_c = 25.0
        RH = 50.0

        d_min = renaudo_extinction_diameter(H, T_c, RH)
        # For agricultural sprays at 50 cm, d_min is typically in 30 - 150 um range
        @test d_min > 30e-6
        @test d_min < 160e-6

        # Higher release height allows more flight time -> larger droplets vanish completely
        d_min_high = renaudo_extinction_diameter(1.0, T_c, RH)
        @test d_min_high > d_min
        # Ratio scales as (H2/H1)^(1/4)
        @test isapprox(d_min_high / d_min, 2.0^0.25; atol=1e-3)
    end

    # 3. Renaudo Deposition Time & Landing Diameter
    @testset "Renaudo In-Flight Deposition & Landing Diameter" begin
        H = 0.5
        T_c = 25.0
        RH = 50.0
        d_min = renaudo_extinction_diameter(H, T_c, RH)

        # Droplet smaller than d_min evaporates completely in flight
        d_sub = 0.8 * d_min
        @test isinf(renaudo_deposition_time(d_sub, H, T_c, RH))
        @test renaudo_landing_diameter(d_sub, H, T_c, RH) == 0.0

        # Droplet larger than d_min lands with remaining mass
        d_super = 150e-6
        t_dep = renaudo_deposition_time(d_super, H, T_c, RH)
        d_land = renaudo_landing_diameter(d_super, H, T_c, RH)

        @test isfinite(t_dep)
        @test t_dep > 0.0
        @test d_land > 0.0
        @test d_land < d_super

        # Very large droplet (e.g. 500 um) lands with almost full size
        d_huge = 500e-6
        d_land_huge = renaudo_landing_diameter(d_huge, H, T_c, RH)
        @test isapprox(d_land_huge, d_huge; rtol=0.05)
    end

    # 4. Chiu's Group Evaporation Number (Sirignano 2010)
    @testset "Chiu Group Evaporation Number & Regimes" begin
        # Isolated droplet limit
        G_iso = chiu_group_evaporation_number(10, 10e-6, 0.1) # 10 droplets in 10 cm cloud
        @test G_iso < 0.1
        @test chiu_evaporation_regime(G_iso) == :single_droplet

        # Intermediate internal group regime
        G_inter = chiu_group_evaporation_number(10_000, 20e-6, 0.05)
        @test G_inter >= 0.1
        @test G_inter <= 100.0
        @test chiu_evaporation_regime(G_inter) == :internal_group

        # Dense spray sheath regime (G > 100)
        G_dense = chiu_group_evaporation_number(2_000_000, 30e-6, 0.01)
        @test G_dense > 100.0
        @test chiu_evaporation_regime(G_dense) == :external_group
    end
end
