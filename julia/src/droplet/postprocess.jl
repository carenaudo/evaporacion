"""
    postprocess.jl

Vectorized post-processing of droplet trajectories into typed DataFrames.
"""

"""
    postprocess_trajectory(sol, params::DropletParams) -> DataFrame

Evaluates all diagnostic algebraic variables at every saved time point in `sol`.
"""
function postprocess_trajectory(sol, params::DropletParams)
    n = length(sol.t)
    df = DataFrame(
        t = sol.t,
        m = [sol.u[i][1] for i in 1:n],
        T_bar = [sol.u[i][2] for i in 1:n],
        Ts = zeros(n),
        d = zeros(n),
        r_s = zeros(n),
        T_eg = zeros(n),
        m_dot = zeros(n),
        Re = zeros(n),
        Pr = zeros(n),
        Sc = zeros(n),
        Nu0 = zeros(n),
        Sh0 = zeros(n),
        Nu_star = zeros(n),
        Sh_star = zeros(n),
        B_M = zeros(n),
        B_T = zeros(n),
        F_M = zeros(n),
        F_T = zeros(n),
        chi_l = zeros(n),
        chi_m = zeros(n),
        U_s = zeros(n),
        h = zeros(n),
        q_g = zeros(n),
        q_l = zeros(n)
    )

    for i in 1:n
        m_i = sol.u[i][1]
        T_bar_i = sol.u[i][2]
        if m_i > 1e-18
            # If solution is from Strategy C (:dae), Ts is already directly available in sol.u[i][3]
            Ts_i = (params.ts_solver === :dae && length(sol.u[i]) >= 3) ? sol.u[i][3] : solve_surface_temperature(m_i, T_bar_i, params)
            alg = calc_droplet_algebraic(m_i, T_bar_i, Ts_i, params)
            df.Ts[i] = Ts_i
            df.d[i] = alg.d
            df.r_s[i] = alg.r_s
            df.T_eg[i] = alg.T_eg
            df.m_dot[i] = alg.m_dot
            df.Re[i] = alg.Re
            df.Pr[i] = alg.Pr
            df.Sc[i] = alg.Sc
            df.Nu0[i] = alg.Nu0
            df.Sh0[i] = alg.Sh0
            df.Nu_star[i] = alg.Nu_star
            df.Sh_star[i] = alg.Sh_star
            df.B_M[i] = alg.B_M
            df.B_T[i] = alg.B_T
            df.F_M[i] = alg.F_M
            df.F_T[i] = alg.F_T
            df.chi_l[i] = alg.chi_l
            df.chi_m[i] = alg.chi_m
            df.U_s[i] = alg.U_s
            df.h[i] = alg.h
            df.q_g[i] = alg.q_g
            df.q_l[i] = alg.q_l
        end
    end

    return df
end
