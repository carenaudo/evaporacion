"""
    ensemble.jl

Multithreaded parametric sweep runner using SciML EnsembleProblem.
"""

"""
    run_parametric_sweep(configs::Vector{SurrogateConfig}; solver=Tsit5(), batch_size=100) -> DataFrame

Executes a parametric sweep of surrogate runs in parallel across all available Julia CPU threads.
"""
function run_parametric_sweep(configs::Vector{SurrogateConfig}; solver=Tsit5(), batch_size=100)
    N = length(configs)
    if N == 0
        return DataFrame()
    end

    cfg0 = configs[1]
    V0 = (π / 6.0) * (cfg0.d0_mm^3)
    u0 = cfg0.model_type in (:modelo_8, :modelo_8_dynamic) ? [V0, cfg0.T0] : [V0]
    cb = ContinuousCallback((u, t, integrator) -> u[1] - 1e-9, terminate!)

    prob0 = ODEProblem(surrogate_ode_rhs!, u0, (0.0, cfg0.t_end), cfg0)

    prob_func = function (prob, i, repeat)
        cfg_i = configs[i]
        V0_i = (π / 6.0) * (cfg_i.d0_mm^3)
        u0_i = cfg_i.model_type in (:modelo_8, :modelo_8_dynamic) ? [V0_i, cfg_i.T0] : [V0_i]
        remake(prob, u0=u0_i, tspan=(0.0, cfg_i.t_end), p=cfg_i)
    end

    ensemble_prob = EnsembleProblem(prob0; prob_func=prob_func)
    sim = solve(ensemble_prob, solver, EnsembleThreads(); trajectories=N, batch_size=batch_size, callback=cb)

    # Collect metrics into a DataFrame
    df = DataFrame(
        case_id = 1:N,
        model_type = [String(c.model_type) for c in configs],
        d0_mm = [c.d0_mm for c in configs],
        T_inf = [c.T_inf for c in configs],
        T_wb = [c.T_wb for c in configs],
        U = [c.U for c in configs],
        HR = [c.HR for c in configs],
        t_evap_s = zeros(N),
        final_V_mm3 = zeros(N)
    )

    for i in 1:N
        sol_i = sim[i]
        df.t_evap_s[i] = sol_i.t[end]
        df.final_V_mm3[i] = sol_i.u[end][1]
    end

    return df
end
