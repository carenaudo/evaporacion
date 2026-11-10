using Test
using Evaporacion
using Evaporacion.Droplet
using OrdinaryDiffEq
using DataFrames

@testset "Surface Temperature (T_s) Solver Strategies & Parity" begin
    T_inf = 25.0
    Y_inf = 0.01
    P_atm = 101325.0
    d0 = 200e-6
    rho_l0 = Evaporacion.Droplet.liquid_density_water(T_inf + 273.15)
    m0 = rho_l0 * (π / 6.0) * d0^3
    T0 = 20.0

    # 1. Pointwise surface temperature evaluation parity
    p_brent = DropletParams(T_inf=T_inf, Y_inf=Y_inf, P_atm=P_atm, profile=:parabolico, ts_solver=:brent)
    p_newton = DropletParams(T_inf=T_inf, Y_inf=Y_inf, P_atm=P_atm, profile=:parabolico, ts_solver=:newton)
    p_halley = DropletParams(T_inf=T_inf, Y_inf=Y_inf, P_atm=P_atm, profile=:parabolico, ts_solver=:halley)
    p_biot = DropletParams(T_inf=T_inf, Y_inf=Y_inf, P_atm=P_atm, profile=:parabolico, ts_solver=:analytical_biot)
    p_halley_ana = DropletParams(T_inf=T_inf, Y_inf=Y_inf, P_atm=P_atm, profile=:parabolico, ts_solver=:analytical_halley)
    p_cheb = DropletParams(T_inf=T_inf, Y_inf=Y_inf, P_atm=P_atm, profile=:parabolico, ts_solver=:chebyshev)
    p_lut = DropletParams(T_inf=T_inf, Y_inf=Y_inf, P_atm=P_atm, profile=:parabolico, ts_solver=:lut)

    Ts_brent = solve_surface_temperature(m0, T0, p_brent)
    Ts_newton = solve_surface_temperature(m0, T0, p_newton)
    Ts_halley = solve_surface_temperature(m0, T0, p_halley)
    Ts_biot = solve_surface_temperature(m0, T0, p_biot)
    Ts_halley_ana = solve_surface_temperature(m0, T0, p_halley_ana)
    Ts_cheb = solve_surface_temperature(m0, T0, p_cheb)
    Ts_lut = solve_surface_temperature(m0, T0, p_lut)

    # Sanity bounds
    @test 10.0 < Ts_brent < T_inf

    # Strict point parity (in Kelvin / Celsius)
    @test isapprox(Ts_newton, Ts_brent, atol=1e-5)      # < 0.01 mK
    @test isapprox(Ts_halley, Ts_brent, atol=1e-5)      # < 0.01 mK
    @test isapprox(Ts_halley_ana, Ts_brent, atol=1e-4)  # < 0.1 mK
    @test isapprox(Ts_biot, Ts_brent, atol=1e-3)        # < 1.0 mK
    @test isapprox(Ts_cheb, Ts_brent, atol=1e-3)        # < 1.0 mK
    @test isapprox(Ts_lut, Ts_brent, atol=1e-3)         # < 1.0 mK

    # 2. Invariant Liquid State Precomputation
    liq = precompute_liquid_state(m0, T0, p_newton)
    @test liq.d > 0.0
    @test liq.r_s == liq.d / 2.0
    @test liq.rho_l > 900.0
    @test liq.cp_l > 4000.0
    @test liq.k_l > 0.5

    # 3. Quiescent 1D Collapse & Surrogate Structure
    p_quiescent = DropletParams(T_inf=T_inf, Y_inf=Y_inf, P_atm=P_atm, U_rel=0.0, terminal_velocity=false)
    surr_1d = build_chebyshev_ts_surrogate(p_quiescent)
    @test surr_1d.is_1d == true
    @test size(surr_1d.coeffs, 1) == 1

    p_falling = DropletParams(T_inf=T_inf, Y_inf=Y_inf, P_atm=P_atm, terminal_velocity=true)
    surr_2d = build_chebyshev_ts_surrogate(p_falling)
    @test surr_2d.is_1d == false
    @test size(surr_2d.coeffs, 1) > 1

    # 4. Trajectory Evaporation Time Parity across ODE and DAE Solvers
    tspan = (0.0, 50.0)
    sol_brent = solve_droplet(p_brent, [m0, T0], tspan)
    sol_newton = solve_droplet(p_newton, [m0, T0], tspan)
    sol_halley = solve_droplet(p_halley, [m0, T0], tspan)
    sol_biot = solve_droplet(p_biot, [m0, T0], tspan)
    sol_halley_ana = solve_droplet(p_halley_ana, [m0, T0], tspan)
    sol_cheb = solve_droplet(p_cheb, [m0, T0], tspan)
    sol_lut = solve_droplet(p_lut, [m0, T0], tspan)

    p_dae = DropletParams(T_inf=T_inf, Y_inf=Y_inf, P_atm=P_atm, profile=:parabolico, ts_solver=:dae)
    sol_dae = solve_droplet(p_dae, [m0, T0], tspan)

    t_evap_brent = sol_brent.t[end]
    @test t_evap_brent > 1.0

    # All methods agree on evaporation lifetime within 0.05%
    @test isapprox(sol_newton.t[end], t_evap_brent, rtol=1e-4)
    @test isapprox(sol_halley.t[end], t_evap_brent, rtol=1e-4)
    @test isapprox(sol_biot.t[end], t_evap_brent, rtol=1e-4)
    @test isapprox(sol_halley_ana.t[end], t_evap_brent, rtol=1e-4)
    @test isapprox(sol_cheb.t[end], t_evap_brent, rtol=1e-4)
    @test isapprox(sol_lut.t[end], t_evap_brent, rtol=1e-4)
    @test isapprox(sol_dae.t[end], t_evap_brent, rtol=5e-4)

    # 5. Dual Sensitivity Propagation
    import ForwardDiff
    dual_m = ForwardDiff.Dual(m0, 1.0)
    Ts_dual = solve_surface_temperature(dual_m, T0, p_cheb)
    @test Ts_dual isa ForwardDiff.Dual
    @test isapprox(ForwardDiff.value(Ts_dual), Ts_cheb, atol=1e-5)
    @test isfinite(ForwardDiff.partials(Ts_dual)[1])

    # 6. Post-processing with DAE trajectory
    df_dae = postprocess_trajectory(sol_dae, p_dae)
    @test nrow(df_dae) == length(sol_dae.t)
    @test df_dae.m_dot[1] > 0.0
    @test all(df_dae.Ts .> 0.0)
end
