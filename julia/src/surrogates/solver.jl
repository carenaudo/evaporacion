"""
    solver.jl

ODE integrator for all empirical and physics-informed surrogate models.
"""

"""
    SurrogateConfig

Configuration for empirical volume models.
"""
Base.@kwdef struct SurrogateConfig
    model_type::Symbol = :modelo_6b_refit
    d0_mm::Float64 = 1.0
    T_inf::Float64 = 25.0
    T_wb::Float64 = 18.0
    U::Float64 = 0.1
    HR::Float64 = 0.5
    T0::Float64 = 22.0
    t_end::Float64 = 30000.0
    dynamic_params::Union{Nothing, Vector{Float64}} = nothing
    coeffs::Any = nothing
end

"""
    surrogate_ode_rhs!(du, u, p::SurrogateConfig, t)

In-place ODE RHS for 1-equation and 2-equation surrogate models.
"""
function surrogate_ode_rhs!(du, u, p::SurrogateConfig, t)
    V = u[1] # Volume in mm^3
    if V <= 1e-9
        du[1] = 0.0
        return nothing
    end

    # 1-equation models
    d_mm = sphere_diameter_from_volume(V)
    C = if p.model_type === :modelo_4
        calc_C_model4(d_mm, p.U, p.T_inf)
    elseif p.model_type === :modelo_5
        calc_C_model5(d_mm, p.U, p.T_inf)
    elseif p.model_type === :modelo_6
        calc_C_model6(d_mm, p.U, p.T_inf)
    elseif p.model_type === :modelo_6b
        calc_C_model6b(d_mm, p.U, p.T_inf)
    elseif p.model_type === :modelo_6b_refit
        calc_C_model6b_refit(d_mm, p.U, p.T_inf; coeffs=p.coeffs)
    elseif p.model_type === :modelo_7
        calc_C_model7(d_mm, p.U, p.T_inf, p.HR, p.T_wb, p.T0; coeffs=p.coeffs)
    elseif p.model_type === :modelo_2b
        calc_C_model2b(d_mm, p.U, p.T_inf, p.T_wb)
    else
        1.8e-4
    end

    S_mm2 = sphere_surface_from_d(d_mm)
    dT = p.T_inf - p.T_wb
    du[1] = -C * sqrt(S_mm2) * dT

    return nothing
end

"""
    solve_surrogate(cfg::SurrogateConfig; alg=Tsit5(), rtol=1e-6, atol=1e-9)

Solve the surrogate evaporation trajectory. Returns an ODESolution.
"""
function solve_surrogate(cfg::SurrogateConfig; solver=Tsit5(), rtol=1e-6, atol=1e-9)
    V0 = (π / 6.0) * (cfg.d0_mm^3)
    u0 = [V0]

    # Continuous termination callback when volume reaches zero threshold
    cb = ContinuousCallback((u, t, integrator) -> u[1] - 1e-9, terminate!)

    prob = ODEProblem(surrogate_ode_rhs!, u0, (0.0, cfg.t_end), cfg)
    return solve(prob, solver; callback=cb, reltol=rtol, abstol=atol)
end
