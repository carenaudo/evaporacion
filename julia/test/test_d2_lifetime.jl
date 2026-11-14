# An isothermal D² law has constant K, so d^2 = d0^2 - K t exactly. The ODE solution must
# reproduce it for small drops too, whose mass (6.5e-11 kg at 50 µm) is far below any
# absolute tolerance expressed in kg. The solver stops when the mass reaches 1e-15 kg, so
# the reference is the time to that mass, not d0^2/K (they differ by 6e-4 at 50 µm).
@testset "Isothermal D² lifetime matches the analytic solution" begin
    T_K, HR = 283.15, 0.20
    T_inf = T_K - 273.15
    pv = HR * psat_water(T_K)
    M_w, M_a = Evaporacion.Common.Constants.M_WATER, Evaporacion.Common.Constants.M_AIR
    Y_inf = pv * M_w / (pv * M_w + (101325.0 - pv) * M_a)
    T_wb = stull_twb(T_inf, HR * 100)

    for variant in (:d2_classic, :d2_stefan), d0 in (50e-6, 200e-6, 1200e-6)
        p = DropletParams(T_inf = T_inf, Y_inf = Y_inf, model_variant = variant,
                          isothermal = true, T_d_fixed = T_wb)
        m0 = liquid_density_water(T_wb + 273.15) * π / 6 * d0^3
        alg = calc_droplet_algebraic(m0, T_wb, T_wb, p)
        K = 4 * alg.m_dot / (π * alg.rho_l * alg.d)
        d_stop = (6 * 1e-15 / (π * alg.rho_l))^(1 / 3)
        t_exact = (d0^2 - d_stop^2) / K

        sol = solve_droplet(p, [m0, T_wb], (0.0, 3 * t_exact))
        @test string(sol.retcode) == "Terminated"     # stopped by the mass callback, not by tspan
        @test isapprox(sol.t[end], t_exact; rtol = 1e-4)
    end
end
