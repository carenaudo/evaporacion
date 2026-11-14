using Test
using Evaporacion
using Evaporacion.Common
using Evaporacion.Droplet
using Evaporacion.Surrogates
using Evaporacion.Sessile
using Evaporacion.RadialPDE
using DataFrames
using Random
using Lux

@testset "Evaporacion.jl Test Suite" begin

    @testset "Common Properties" begin
        # Saturation vapor pressure at 25 °C (298.15 K)
        p_sat_25 = psat_water(298.15)
        @test isapprox(p_sat_25, 3168.53, rtol=1e-3)

        # Density of liquid water at 25 °C
        rho_25 = liquid_density_water(298.15)
        @test isapprox(rho_25, 997.048, atol=0.01)

        # Latent heat of vaporization at boiling point (373.15 K)
        L_boil = latent_heat_water(373.15)
        @test isapprox(L_boil, 2.257e6, rtol=1e-3)

        # Stull wet-bulb temperature
        twb = stull_twb(35.0, 40.0)
        @test isapprox(twb, 24.514, atol=0.05)
    end

    @testset "Common Correlations & Geometry" begin
        # Blowing factor series limit at B=0
        @test isapprox(blowing_factor(0.0), 1.0, atol=1e-6)

        # Exact values of Wilson & D'Ambrosio g(θ)
        @test isapprox(g_theta_fast(0.5 * π), 2.0, atol=1e-2)
        @test isapprox(g_theta_exact(0.5 * π), 2.0, atol=1e-4)

        # Volume derivatives
        R = 0.001
        th = 0.5 * π
        dV_dth = dV_dtheta(R, th)
        @test dV_dth > 0.0
    end

    @testset "Rigorous Spherical Droplet Model" begin
        params = DropletParams(
            T_inf = 35.0,
            U_rel = 1.0,
            Y_inf = 0.01,
            profile = :parabolico
        )

        m0 = 997.0 * (π / 6.0) * (0.001^3)
        T0 = 22.0

        # Surface temperature root-finding
        Ts = solve_surface_temperature(m0, T0, params)
        @test 10.0 < Ts < 35.0

        alg = calc_droplet_algebraic(m0, T0, Ts, params)
        @test alg.m_dot > 0.0
        @test alg.Re > 0.0
        @test alg.Nu_star >= 2.0
        @test alg.Sh_star >= 2.0

        # Short integration
        sol = solve_droplet(params, [m0, T0], (0.0, 2.0))
        @test length(sol.t) > 2
        @test sol.u[end][1] < sol.u[1][1] # mass decreased

        df = postprocess_trajectory(sol, params)
        @test nrow(df) == length(sol.t)
        @test df.m_dot[1] > 0.0
    end

    @testset "Calibrated Surrogates" begin
        cfg_6b = SurrogateConfig(
            model_type = :modelo_6b_refit,
            d0_mm = 1.0,
            T_inf = 35.0,
            T_wb = 23.0,
            U = 0.5,
            HR = 0.4,
            t_end = 5.0
        )
        sol_6b = solve_surrogate(cfg_6b)
        @test sol_6b.u[end][1] < sol_6b.u[1][1]

        cfg_7 = SurrogateConfig(
            model_type = :modelo_7,
            d0_mm = 1.0,
            T_inf = 35.0,
            T_wb = 23.0,
            U = 0.5,
            HR = 0.4,
            T0 = 22.0,
            t_end = 5.0
        )
        sol_7 = solve_surrogate(cfg_7)
        @test length(sol_7.u[1]) == 1
        @test sol_7.u[end][1] < sol_7.u[1][1]
    end

    @testset "Sessile Droplet Suite" begin
        cfg_cr = SessileDropletConfig(
            R0 = 1e-3,
            theta0 = 0.5 * π,
            mode = :CR,
            t_end = 10.0
        )
        t_cr_ana = lifetime_cr(cfg_cr)
        @test t_cr_ana > 0.0

        sol_cr = solve_sessile(cfg_cr)
        @test sol_cr.theta[end] < sol_cr.theta[1]

        cfg_ca = SessileDropletConfig(
            R0 = 1e-3,
            theta0 = 0.5 * π,
            mode = :CA,
            t_end = 10.0
        )
        t_ca_ana = lifetime_ca(cfg_ca)
        @test t_ca_ana > 0.0

        sol_ca = solve_sessile(cfg_ca)
        @test sol_ca.R[end] < sol_ca.R[1]
    end

    @testset "1D Radial Moving-Boundary PDE" begin
        cfg_pde = RadialPDEConfig(
            n_nodes = 11,
            d0 = 0.001,
            T_drop_0 = 22.0,
            T_air = 35.0,
            RH = 0.4,
            t_end = 2.0
        )
        sol_pde = solve_radial_pde(cfg_pde)
        @test sol_pde.R[end] < sol_pde.R[1]
        @test size(sol_pde.T_profile, 2) == 11
    end

    @testset "SciML Neural ODE Models" begin
        import Random
        rng = Random.default_rng()

        # Option 3b MLP Model
        model_3b = Option3bNeuralODE()
        ps_3b, st_3b = Lux.setup(rng, model_3b)
        u0_norm = Float32[0.0, 0.0]
        cond_norm = Float32[0.0, 0.0, 0.0, 0.0, 0.0]
        tsteps = Float32[0.0, 0.5, 1.0]

        sol_3b = solve_neural_ode_3b(model_3b, ps_3b, st_3b, u0_norm, cond_norm, tsteps)
        @test length(sol_3b.u) == 3

        # Option 3c Hybrid Attention Model
        model_3c = Option3cHybridODE()
        ps_3c, st_3c = Lux.setup(rng, model_3c)
        state_cond = zeros(Float32, 7)
        delta_u = Evaporacion.Neural.Option3cHybridAttn.evaluate_hybrid_residual(model_3c, ps_3c, st_3c, state_cond)
        @test length(delta_u) == 2
    end

    @testset "Physics-Informed Neural & Operator Surrogates" begin
        rng = Random.MersenneTwister(123)
        p_c    = PINNCSurrogate(PINNSurrogates.init_mlp_params(rng, [9, 24, 24, 1]; out_scale=0.01))
        p_2eq  = PINN2EqSurrogate(PINNSurrogates.init_mlp_params(rng, [10, 32, 32, 3]; out_scale=0.01))
        p_deep = PIDeepONetSurrogate(
            PINNSurrogates.init_mlp_params(rng, [10, 48, 48, 3 + 3 * 16]; out_scale=0.02),
            PINNSurrogates.init_mlp_params(rng, [4, 24, 24, 16]; out_scale=0.05)
        )
        p_res  = ResMLPSurrogate(
            PINNSurrogates.init_mlp_params(rng, [5, 64, 64, 64, 1]; out_scale=0.05),
            PINNSurrogates.init_mlp_params(rng, [5, 48, 48, 2]; out_scale=0.05),
            zeros(5), ones(5)
        )

        # 1. PINN-C ODE vs Gauss-Legendre consistency & monotonic evaporation
        t_gl  = PINNSurrogates.pinn_c_predict_tevap(p_c, 0.5, 1.0, 15.0, 35.0, 0.4)
        sol_c = solve_pinn_c(p_c, 0.5, 1.0, 15.0, 35.0, 0.4; t_end=30000.0)
        @test t_gl > 0.0
        @test isapprox(Float64(sol_c.t[end]), t_gl; rtol=1e-2)
        @test sol_c.u[end][1] < sol_c.u[1][1]

        # 2. PINN-2Eq coupled [V, T_bar] conservation
        sol_2eq = solve_pinn_2eq(p_2eq, 0.5, 1.0, 15.0, 35.0, 0.4; t_end=30000.0)
        @test length(sol_2eq.u[1]) == 2
        @test sol_2eq.u[end][1] < sol_2eq.u[1][1]
        @test 10.0 < sol_2eq.u[end][2] < 35.0

        # 3. PI-DeepONet exact Hard Boundary Conditions: d(0)=d0, d(1)=0, T(0)=T0
        # For T0=15°C < T_dew(~19.4°C at Tinf=35°C, HR=40%), initial condensation occurs (m_dot_0 < 0);
        # for T0=25°C > T_dew, initial evaporation occurs (m_dot_0 > 0).
        op = deeponet_eval_case(p_deep, 0.5, 1.0, 15.0, 35.0, 0.4; taus=[0.0, 0.5, 1.0])
        @test isapprox(op.d_mm[1], 0.5; atol=1e-12)
        @test isapprox(op.d_mm[end], 0.0; atol=1e-12)
        @test isapprox(op.T_bar[1], 15.0; atol=1e-12)
        @test op.m_dot_0_kgs < 0.0
        @test op.t_evap_s > 0.0
        op_evap = deeponet_eval_case(p_deep, 0.5, 1.0, 25.0, 35.0, 0.4; taus=[0.0, 0.5, 1.0])
        @test op_evap.m_dot_0_kgs > 0.0

        # 4. Bundle serialization round-trip
        tmp_npz = tempname() * ".npz"
        bundle = PINNSurrogateBundle(p_c, p_2eq, p_deep, p_res)
        save_pinn_bundle(bundle, tmp_npz)
        loaded = load_pinn_bundle(tmp_npz)
        op2 = deeponet_eval_case(loaded.pi_deeponet, 0.5, 1.0, 15.0, 35.0, 0.4; taus=[0.0, 0.5, 1.0])
        @test isapprox(op2.t_evap_s, op.t_evap_s; rtol=1e-12)
        rm(tmp_npz; force=true)
    end

    include("test_phase1_evaporation.jl")
    include("test_phase2_sessile.jl")
    include("test_phase3_kinetic.jl")
    include("test_phase4_multicomponent_saline.jl")
    include("test_phase5_extinction.jl")
    include("test_audit_fixes.jl")
    include("test_ts_solvers.jl")
    include("test_d2_lifetime.jl")

end
