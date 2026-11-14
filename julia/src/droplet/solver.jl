"""
    solver.jl

Differential equation solver interface for the rigorous spherical droplet model.
Uses OrdinaryDiffEq stiff/non-stiff algorithms and continuous callbacks.
"""

"""
    droplet_ode_rhs!(du, u, p::DropletParams, t)

ODE right-hand side for the spherical droplet:
    u[1] = m      (droplet mass [kg])
    u[2] = T_bar  (droplet mean temperature [°C])
"""
function droplet_ode_rhs!(du, u, p::DropletParams, t)
    m     = u[1]
    T_bar = u[2]

    if m <= 1e-18
        du[1] = zero(m)
        du[2] = zero(T_bar)
        return nothing
    end

    if p.isothermal
        Ts = isnan(p.T_d_fixed) ? T_bar : p.T_d_fixed
        alg = calc_droplet_algebraic(m, Ts, Ts, p)
        du[1] = -alg.m_dot
        du[2] = zero(T_bar)
        return nothing
    end

    # Implicit surface temperature and algebraic diagnostics in a single unified pass
    Ts, alg = solve_droplet_state(m, T_bar, p)

    # dm/dt = -m_dot
    du[1] = -alg.m_dot

    # dT_bar/dt with Zhou et al. (2017) moving-boundary term and Angelilli et al. (2024) variable-cp factor
    term1 = 3.0 * alg.h * (p.T_inf - Ts) / (alg.rho_l * alg.cp_l * alg.r_s)
    term2 = 3.0 * alg.r_dot_s * alg.DHev / (alg.r_s * alg.cp_l)
    term3 = p.zhou_moving_boundary ? (3.0 * alg.r_dot_s * (Ts - T_bar) / alg.r_s) : zero(term1)
    
    dT_dt_raw = term1 + term2 + term3

    if p.variable_cp_factor
        t_C = T_bar
        dcp_dT = -2.0 + 0.02 * t_C
        f_C = 1.0 + ((T_bar + 273.15) / alg.cp_l) * dcp_dT
        du[2] = dT_dt_raw / f_C
    else
        du[2] = dT_dt_raw
    end

    return nothing
end

"""
    droplet_dae_rhs!(du, u, p::DropletParams, t)

Index-1 Differential-Algebraic Equation (DAE) right-hand side with singular mass matrix:
    M * du/dt = f(u, p, t)
where M = diag(1, 1, 0) and state vector:
    u[1] = m      (differential: dm/dt = -m_dot)
    u[2] = T_bar  (differential: dT_bar/dt = ...)
    u[3] = Ts     (algebraic constraint: 0 = intensive_surface_residual(Ts, m, T_bar))
"""
function droplet_dae_rhs!(du, u, p::DropletParams, t)
    m     = u[1]
    T_bar = u[2]
    Ts    = u[3]

    if m <= 1e-18
        du[1] = zero(m)
        du[2] = zero(T_bar)
        du[3] = zero(Ts)
        return nothing
    end

    alg = calc_droplet_algebraic(m, T_bar, Ts, p)

    # 1. dm/dt = -m_dot
    du[1] = -alg.m_dot

    # 2. dT_bar/dt
    term1 = 3.0 * alg.h * (p.T_inf - Ts) / (alg.rho_l * alg.cp_l * alg.r_s)
    term2 = 3.0 * alg.r_dot_s * alg.DHev / (alg.r_s * alg.cp_l)
    term3 = p.zhou_moving_boundary ? (3.0 * alg.r_dot_s * (Ts - T_bar) / alg.r_s) : zero(term1)
    dT_dt_raw = term1 + term2 + term3

    if p.variable_cp_factor
        t_C = T_bar
        dcp_dT = -2.0 + 0.02 * t_C
        f_C = 1.0 + ((T_bar + 273.15) / alg.cp_l) * dcp_dT
        du[2] = dT_dt_raw / f_C
    else
        du[2] = dT_dt_raw
    end

    # 3. Algebraic surface energy balance residual (intensive form normalized by 4π*r_s)
    # R_tilde(Ts) = [q_g - m_dot * DHev - q_l] / (4π * r_s) = 0
    du[3] = (alg.q_g - alg.m_dot * alg.DHev - alg.q_l) / (4.0 * π * alg.r_s)

    return nothing
end

"""
    make_mass_zero_callback()

Continuous callback to terminate integration when droplet mass drops below 1e-15 kg.
"""
function make_mass_zero_callback()
    condition(u, t, integrator) = u[1] - 1e-15
    affect!(integrator) = terminate!(integrator)
    return ContinuousCallback(condition, affect!)
end

"""
    solve_droplet(params::DropletParams, u0, tspan; solver=nothing, rtol=1e-6, atol=1e-9, kwargs...)

Integrate the droplet evaporation ODE/DAE system from tspan[1] to tspan[2].
Automatically handles Strategy C (:dae) with singular mass-matrix formulation.

A scalar `atol` is applied to the mass relative to the initial mass (`atol * m0`, in kg)
and to `T_bar` as is (°C). Droplet masses are 1e-13 to 1e-6 kg, so an absolute tolerance in
kg would exceed the whole mass of a small drop and leave its error uncontrolled. Pass a
vector to set both tolerances explicitly.
"""
function solve_droplet(
    params::DropletParams, u0::AbstractVector, tspan::Tuple{Real, Real};
    solver=nothing, rtol=1e-6, atol=1e-9, kwargs...
)
    cb = make_mass_zero_callback()

    if params.ts_solver === :dae
        # Strategy C: Index-1 DAE with mass matrix M = diag(1, 1, 0)
        M = [1.0 0.0 0.0;
             0.0 1.0 0.0;
             0.0 0.0 0.0]
        
        # Consistent initial Ts
        Ts0 = length(u0) >= 3 ? Float64(u0[3]) : Float64(solve_surface_temperature(u0[1], u0[2], params))
        u0_dae = [Float64(u0[1]), Float64(u0[2]), Ts0]

        dae_func = ODEFunction(droplet_dae_rhs!; mass_matrix=M)
        prob = ODEProblem(dae_func, u0_dae, float.(tspan), params)

        # Default stiff DAE solver if not specified (autodiff=false avoids Dual instabilities).
        # Use a scalar abstol so NonlinearSolveBase DAE initialization does not error on vector tolerances.
        dae_solver = solver === nothing ? Rodas5P(autodiff=false) : solver
        dae_atol = atol isa Real ? min(Float64(atol) * 1e-4, 1e-13) : atol
        return solve(prob, dae_solver; callback=cb, reltol=rtol, abstol=dae_atol, kwargs...)
    else
        # Standard 2-variable ODE [m, T_bar]
        prob = ODEProblem(droplet_ode_rhs!, float.(u0[1:2]), float.(tspan), params)
        ode_solver = solver === nothing ? Tsit5() : solver
        abstol = atol isa Real ? [Float64(atol) * abs(Float64(ForwardDiff.value(u0[1]))), Float64(atol)] : atol
        return solve(prob, ode_solver; callback=cb, reltol=rtol, abstol=abstol, kwargs...)
    end
end
