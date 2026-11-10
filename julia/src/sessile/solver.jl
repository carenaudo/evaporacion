"""
    solver.jl

Solvers and exact analytical lifetimes for sessile droplet evaporation modes
(CR, CA, SS, and SJ) based on Wilson & D'Ambrosio (2023) Tables 1 & 2.
"""

"""
    SessileDropletConfig

Configuration for a sessile droplet evaporation simulation.
"""
Base.@kwdef struct SessileDropletConfig
    D::Float64 = 2.6e-5            # Binary diffusivity [m^2/s]
    delta_c::Float64 = 0.02        # (c_sat - c_inf) [kg/m^3]
    rho_l::Float64 = 997.0         # Liquid density [kg/m^3]
    R0::Float64 = 1.0e-3           # Initial contact radius [m]
    theta0::Float64 = 0.5 * π      # Initial contact angle [rad]
    theta_rec::Float64 = 0.1 * π   # Receding angle for Stick-Slip mode [rad]
    theta_min::Float64 = 0.1 * π   # Min angle for Stick-Jump mode [rad]
    theta_max::Float64 = 0.5 * π   # Max angle for Stick-Jump mode [rad]
    mode::Symbol = :CR             # :CR, :CA, :SS, :SJ
    t_end::Float64 = 1000.0
    thermal::Union{Nothing, SessileThermalParams} = nothing
end

struct SessileSolution
    t::Vector{Float64}
    R::Vector{Float64}
    theta::Vector{Float64}
    V::Vector{Float64}
    mode::Symbol
    t_evap_analytical::Float64
end

# =============================================================================
# EXACT ANALYTICAL LIFETIMES (Wilson & D'Ambrosio 2023 Table 1)
# =============================================================================

"""
    lifetime_cr(cfg::SessileDropletConfig) -> Float64

Exact lifetime in Constant Radius (pinned) mode:
    t_CR = (ρ * R₀^2 / (D * Δc)) * ∫_0^θ₀ (1 / g(θ)) dθ
"""
function lifetime_cr(cfg::SessileDropletConfig)
    integral_val, _ = quadgk(th -> 1.0 / g_theta_fast(th), 1e-6, cfg.theta0; rtol=1e-6)
    return (cfg.rho_l * (cfg.R0^2) / (cfg.D * cfg.delta_c)) * integral_val
end

"""
    lifetime_ca(cfg::SessileDropletConfig) -> Float64

Exact lifetime in Constant Angle (unpinned) mode:
    t_CA = (ρ * R₀^2 / (2 * D * Δc)) * (sin(θ₀) * (2 + cos(θ₀)) / g(θ₀))
"""
function lifetime_ca(cfg::SessileDropletConfig)
    s = sin(cfg.theta0)
    c = cos(cfg.theta0)
    g_val = g_theta_fast(cfg.theta0)
    geom_factor = s * (2.0 + c) / g_val
    return (cfg.rho_l * (cfg.R0^2) / (2.0 * cfg.D * cfg.delta_c)) * geom_factor
end

"""
    lifetime_ss(cfg::SessileDropletConfig) -> Float64

Exact lifetime in Stick-Slip (pinned until θ_rec, then unpinned) mode.
"""
function lifetime_ss(cfg::SessileDropletConfig)
    th_star = cfg.theta_rec
    integral_val, _ = quadgk(th -> 2.0 / g_theta_fast(th), th_star, cfg.theta0; rtol=1e-6)
    s_star = sin(th_star)
    c_star = cos(th_star)
    g_star = g_theta_fast(th_star)
    geom_factor = s_star * (2.0 + c_star) / g_star
    return (cfg.rho_l * (cfg.R0^2) / (2.0 * cfg.D * cfg.delta_c)) * (integral_val + geom_factor)
end

"""
    lifetime_sj(cfg::SessileDropletConfig) -> Float64

Exact lifetime in Stick-Jump (periodic depinning) mode.
"""
function lifetime_sj(cfg::SessileDropletConfig)
    th_min = cfg.theta_min
    th_max = cfg.theta_max
    integral_val, _ = quadgk(th -> 1.0 / g_theta_fast(th), th_min, th_max; rtol=1e-6)

    s_min, c_min = sin(th_min), cos(th_min)
    s_max, c_max = sin(th_max), cos(th_max)

    lambda = (s_min / s_max) * ((2.0 + c_min) / (2.0 + c_max)) * (((1.0 + c_max) / (1.0 + c_min))^2)
    series_factor = 1.0 / (1.0 - lambda^(2.0 / 3.0))

    return (cfg.rho_l * (cfg.R0^2) / (cfg.D * cfg.delta_c)) * integral_val * series_factor
end

# =============================================================================
# THIN-DROPLET ANALYTICAL LIFETIMES (Wilson & D'Ambrosio 2023 Table 2)
# =============================================================================

"""
    lifetime_cr_thin(R0::Real, theta0::Real, D::Real, delta_c::Real, rho_l::Real) -> Real

Exact thin-droplet (θ₀ ≪ 1) lifetime in Constant Radius (pinned) mode:
    t_CR,thin = π * ρ * θ₀ * R₀^2 / (16 * D * Δc)
"""
function lifetime_cr_thin(R0::Real, theta0::Real, D::Real, delta_c::Real, rho_l::Real)
    return Float64(π * rho_l * theta0 * (R0^2) / (16.0 * D * delta_c))
end

"""
    lifetime_ca_thin(R0::Real, theta0::Real, D::Real, delta_c::Real, rho_l::Real) -> Real

Exact thin-droplet (θ₀ ≪ 1) lifetime in Constant Angle (unpinned) mode:
    t_CA,thin = 3π * ρ * θ₀ * R₀^2 / (32 * D * Δc) = 1.5 * t_CR,thin
"""
function lifetime_ca_thin(R0::Real, theta0::Real, D::Real, delta_c::Real, rho_l::Real)
    return Float64(3.0 * π * rho_l * theta0 * (R0^2) / (32.0 * D * delta_c))
end

"""
    lifetime_ss_thin(R0::Real, theta0::Real, theta_star::Real, D::Real, delta_c::Real, rho_l::Real) -> Real

Exact thin-droplet (θ₀ ≪ 1) lifetime in Stick-Slide mode (verified from Wilson 2023 Table 2):
    t_SS,thin = π * ρ * (2θ₀ + θ*) * R₀^2 / (32 * D * Δc)
"""
function lifetime_ss_thin(R0::Real, theta0::Real, theta_star::Real, D::Real, delta_c::Real, rho_l::Real)
    return Float64(π * rho_l * (2.0 * theta0 + theta_star) * (R0^2) / (32.0 * D * delta_c))
end

"""
    lifetime_sj_thin(R0::Real, theta0::Real, theta_star::Real, n_jumps::Int, D::Real, delta_c::Real, rho_l::Real) -> Real

Exact thin-droplet (θ₀ ≪ 1) lifetime after n stick-jump cycles to extinction (Wilson 2023 Table 2):
For n = 1 (no jumps), t = t_CR,thin.
For n > 1 with jumps, time of (n-1) stick phases plus final CR evaporation to extinction:
    t_n = t_CR,thin * [ (1 - θ*/θ₀) * (1 - λ^(2(n-1)/3)) / (1 - λ^(2/3)) + λ^(2(n-1)/3) ]
As n -> ∞, t_∞ = t_CR,thin * (1 - θ*/θ₀) / (1 - (θ*/θ₀)^(2/3)).
"""
function lifetime_sj_thin(R0::Real, theta0::Real, theta_star::Real, n_jumps::Int, D::Real, delta_c::Real, rho_l::Real)
    t_cr = lifetime_cr_thin(R0, theta0, D, delta_c, rho_l)
    if theta0 <= theta_star || n_jumps <= 1
        return t_cr
    end
    lambda_2 = (theta_star / theta0)^(2.0 / 3.0)
    pow_n = lambda_2^(n_jumps - 1)
    factor = (1.0 - theta_star / theta0) * (1.0 - pow_n) / (1.0 - lambda_2) + pow_n
    return t_cr * factor
end

# =============================================================================
# STAUBER ET AL. (2015) HYDROPHOBIC SPHERE LIMIT & CRITICAL ANGLE
# =============================================================================

"""
    lifetime_hydrophobic_sphere(V0::Real, D::Real, delta_c::Real, rho_l::Real) -> Real

Exact lifetime t_π of an isolated spherical droplet in the superhydrophobic contact angle limit θ -> π
(Stauber et al. 2015, Eq. 4):
    t_π = (ρ / (2 * D * Δc * ln(2))) * (3*V₀ / (4π))^(2/3)
"""
function lifetime_hydrophobic_sphere(V0::Real, D::Real, delta_c::Real, rho_l::Real)
    R_equiv = (3.0 * V0 / (4.0 * π))^(1.0 / 3.0)
    return Float64((rho_l * (R_equiv^2)) / (2.0 * D * delta_c * log(2.0)))
end

"""
    stauber_critical_angle() -> Float64

The Stauber et al. (2015) non-monotonic critical contact angle θ_crit ≈ 148.2° (2.5865 rad)
where sessile droplet lifetime achieves a local minimum before increasing towards the t_π sphere limit.
"""
function stauber_critical_angle()
    return 148.2 * π / 180.0
end

# =============================================================================
# MULTI-DROPLET PAIR & NONPLANAR WELL LIFETIMES (Wilson 2023 Eqs. 38–40)
# =============================================================================

"""
    lifetime_pair_cr(B::Real, R0::Real, theta0::Real, D::Real, delta_c::Real, rho_l::Real) -> Real

Lifetime of a thin droplet pair with center-to-center separation B in CR mode (Wray et al. 2020):
    t_pair = [1 + (2/π) * arcsin(R₀ / B)] * t_CR,thin
"""
function lifetime_pair_cr(B::Real, R0::Real, theta0::Real, D::Real, delta_c::Real, rho_l::Real)
    t_cr = lifetime_cr_thin(R0, theta0, D, delta_c, rho_l)
    ratio = clamp(R0 / max(B, R0 + 1e-12), 0.0, 1.0)
    return t_cr * (1.0 + (2.0 / π) * asin(ratio))
end

"""
    lifetime_pair_ss(B::Real, R0::Real, theta0::Real, theta_star::Real, D::Real, delta_c::Real, rho_l::Real) -> Real

Lifetime of a thin droplet pair with center-to-center separation B in Stick-Slide mode (Wilson 2023 Eq. 38).
"""
function lifetime_pair_ss(B::Real, R0::Real, theta0::Real, theta_star::Real, D::Real, delta_c::Real, rho_l::Real)
    t_ss = lifetime_ss_thin(R0, theta0, theta_star, D, delta_c, rho_l)
    b_val = max(B, R0 + 1e-12)
    arg_asin = clamp(R0 / b_val, 0.0, 1.0)

    term_num = (2.0 * (2.0 * theta0 + theta_star) * (R0^2) - 3.0 * theta_star * (b_val^2)) * asin(arg_asin) +
               3.0 * theta_star * R0 * sqrt(max(b_val^2 - R0^2, 0.0))
    term_den = π * (2.0 * theta0 + theta_star) * (R0^2)

    return t_ss * (1.0 + term_num / term_den)
end

"""
    lifetime_well_cr(H0::Real, n::Real, R0::Real, theta0::Real, D::Real, delta_c::Real, rho_l::Real) -> Real

Lifetime of a thin pinned droplet inside a shallow nonplanar well z = -H₀ * [1 - (r/R₀)^n]
(D'Ambrosio et al. 2021, Wilson & D'Ambrosio 2023 Eqs. 39–40).
"""
function lifetime_well_cr(H0::Real, n::Real, R0::Real, theta0::Real, D::Real, delta_c::Real, rho_l::Real)
    t_cr = lifetime_cr_thin(R0, theta0, D, delta_c, rho_l)
    n_val = max(Float64(n), 0.1)
    if n_val <= 2.0
        factor = 1.0 + 3.0 * n_val * H0 / ((1.0 + n_val) * theta0 * R0)
    else
        # Wilson & D'Ambrosio (2023) Eq. (40): α(n) -> 0.1369 as n -> ∞ with α(2) = 0
        alpha_n = 0.1369 * (n_val - 2.0) / (n_val - 1.0)
        factor = 1.0 + 2.0 * (1.0 + 8.0 * alpha_n) * H0 / (theta0 * R0)
    end
    return t_cr * factor
end

"""
    lifetime_conduction_limited(R0::Real, theta0::Real, rho_l::Real, L::Real, h_s::Real, k_s::Real, delta_c::Real, dcsat_dT::Real) -> Real

Lifetime under extreme substrate thermal resistance (conduction-limited regime, Wilson 2023 Eq. 30):
    t_cond = (ρ * θ₀ * R₀ * L * h_s * (dc_sat/dT)) / (k_s * Δc)
"""
function lifetime_conduction_limited(
    R0::Real, theta0::Real, rho_l::Real, L::Real, h_s::Real, k_s::Real, delta_c::Real, dcsat_dT::Real
)
    return Float64(rho_l * theta0 * R0 * L * h_s * dcsat_dT / (k_s * delta_c))
end

# =============================================================================
# SAXTON ET AL. (2016) GENERALIZED TANNER'S LAW WITH EVAPORATION
# =============================================================================

"""
    saxton_effective_contact_angle(E_cap::Real, n::Int=1) -> Float64

Evaluates the effective microscopic contact angle Φ₀(E, n) from Saxton et al. (2016, Eq. 3.12):
    Φ₀(E, n) ≈ [1 + Φ_∞(n)^(7/2) * E]^(2/7)
where Φ_∞(1) = 1.750 and Φ_∞(2) = 1.939. Accurate within 2% across all E.
"""
function saxton_effective_contact_angle(E_cap::Real, n::Int=1)
    phi_inf = n == 1 ? 1.750 : 1.939
    e_val = max(Float64(E_cap), 0.0)
    return (1.0 + (phi_inf^3.5) * e_val)^(2.0 / 7.0)
end

"""
    saxton_extinction_scaling(T::Real, Tc::Real, E_cap::Real; n::Int=1) -> Float64

Calculates the exact d^(13/7)-law radius evolution S₀(T) near the extinction singularity
(Saxton et al. 2016, Eq. 3.25):
    S₀(T) ∼ [26 / (5 * Φ_∞(n))]^(7/13) * E^(-2/13) * (Tc - T)^(7/13)
"""
function saxton_extinction_scaling(T::Real, Tc::Real, E_cap::Real; n::Int=1)
    phi_inf = n == 1 ? 1.750 : 1.939
    delta_t = max(Float64(Tc - T), 0.0)
    e_val = max(Float64(E_cap), 1e-12)
    pre = (26.0 / (5.0 * phi_inf))^(7.0 / 13.0)
    return pre * (e_val^(-2.0 / 13.0)) * (delta_t^(7.0 / 13.0))
end

"""
    solve_saxton_spreading_evaporation(V_star::Real, E_cap::Real, t_span::Tuple{Real, Real}; n::Int=1, s0_init::Real=1.0)

Integrates the Saxton et al. (2016, Eq. 3.20) generalized Tanner spreading ODE with evaporation:
    dŝ₀/dt̂ = (1/3) * [ (4*V* / (π * ŝ₀^3))^3 - Φ₀(E / ŝ₀^(1/2), n)^3 ]
Returns an OrdinaryDiffEq ODESolution representing dimensionless contact radius ŝ₀(t̂).
"""
function solve_saxton_spreading_evaporation(
    V_star::Real, E_cap::Real, t_span::Tuple{Real, Real};
    n::Int=1, s0_init::Real=1.0, solver=Tsit5(), rtol=1e-6, atol=1e-8
)
    v_star = Float64(V_star)
    e_cap  = Float64(E_cap)

    rhs_saxton = function (du, u, p, t)
        s0 = max(u[1], 1e-4)
        v_term = (4.0 * v_star / (π * (s0^3)))^3
        phi_0 = saxton_effective_contact_angle(e_cap / sqrt(s0), n)
        du[1] = (1.0 / 3.0) * (v_term - phi_0^3)
    end

    prob = ODEProblem(rhs_saxton, [Float64(s0_init)], float.(t_span))
    return solve(prob, solver; reltol=rtol, abstol=atol)
end

# =============================================================================
# DYNAMIC ODE INTEGRATION
# =============================================================================

"""
    solve_sessile(cfg::SessileDropletConfig; solver=Tsit5(), rtol=1e-8, atol=1e-11) -> SessileSolution

Integrates the sessile droplet evolution in the selected mode (:CR, :CA, :SS, :SJ).
"""
function solve_sessile(cfg::SessileDropletConfig; solver=Tsit5(), rtol=1e-8, atol=1e-11)
    if cfg.mode === :CR
        t_ana = lifetime_cr(cfg)
        t_span = (0.0, min(cfg.t_end, 1.05 * t_ana))

        # In CR mode: dθ/dt = - (D * Δc_eff / (ρ * R0^2)) * g(θ)
        rhs_cr = function (du, u, p, t)
            th = u[1]
            if th <= 1e-6
                du[1] = 0.0
                return nothing
            end
            delta_c_eff = p.thermal === nothing ? p.delta_c : cooling_effective_delta_c(p.R0, th, p.D, p.delta_c, p.thermal)
            du[1] = -(p.D * delta_c_eff / (p.rho_l * (p.R0^2))) * g_theta_fast(th)
            return nothing
        end

        cb = ContinuousCallback((u, t, integrator) -> u[1] - 1e-6, terminate!)
        prob = ODEProblem(rhs_cr, [cfg.theta0], t_span, cfg)
        sol = solve(prob, solver; callback=cb, reltol=rtol, abstol=atol)

        t_vec = sol.t
        theta_vec = [u[1] for u in sol.u]
        R_vec = fill(cfg.R0, length(t_vec))
        V_vec = [sessile_volume_from_r(cfg.R0, th) for th in theta_vec]

        return SessileSolution(t_vec, R_vec, theta_vec, V_vec, :CR, t_ana)

    elseif cfg.mode === :CA
        t_ana = lifetime_ca(cfg)
        t_span = (0.0, min(cfg.t_end, 1.05 * t_ana))

        rhs_ca = function (du, u, p, t)
            R = u[1]
            if R <= 1e-9
                du[1] = 0.0
                return nothing
            end
            delta_c_eff = p.thermal === nothing ? p.delta_c : cooling_effective_delta_c(R, p.theta0, p.D, p.delta_c, p.thermal)
            s = sin(p.theta0)
            c = cos(p.theta0)
            geom = g_theta_fast(p.theta0) / (s * (2.0 + c))
            du[1] = -(p.D * delta_c_eff / (p.rho_l * R)) * geom
            return nothing
        end

        cb = ContinuousCallback((u, t, integrator) -> u[1] - 1e-9, terminate!)
        prob = ODEProblem(rhs_ca, [cfg.R0], t_span, cfg)
        sol = solve(prob, solver; callback=cb, reltol=rtol, abstol=atol)

        t_vec = sol.t
        R_vec = [u[1] for u in sol.u]
        theta_vec = fill(cfg.theta0, length(t_vec))
        V_vec = [sessile_volume_from_r(r, cfg.theta0) for r in R_vec]

        return SessileSolution(t_vec, R_vec, theta_vec, V_vec, :CA, t_ana)

    elseif cfg.mode === :SS
        t_ana = lifetime_ss(cfg)
        # Stage 1: CR mode from theta0 down to theta_rec
        rhs_cr = function (du, u, p, t)
            th = u[1]
            if th <= p.theta_rec
                du[1] = 0.0
                return nothing
            end
            delta_c_eff = p.thermal === nothing ? p.delta_c : cooling_effective_delta_c(p.R0, th, p.D, p.delta_c, p.thermal)
            du[1] = -(p.D * delta_c_eff / (p.rho_l * (p.R0^2))) * g_theta_fast(th)
            return nothing
        end

        cb_cr = ContinuousCallback((u, t, integrator) -> u[1] - cfg.theta_rec, terminate!)
        prob_cr = ODEProblem(rhs_cr, [cfg.theta0], (0.0, cfg.t_end), cfg)
        sol_cr = solve(prob_cr, solver; callback=cb_cr, reltol=rtol, abstol=atol)

        t1 = sol_cr.t[end]
        t_rem = max(0.0, cfg.t_end - t1)

        # Stage 2: CA mode at theta_rec
        rhs_ca = function (du, u, p, t)
            R = u[1]
            if R <= 1e-9
                du[1] = 0.0
                return nothing
            end
            delta_c_eff = p.thermal === nothing ? p.delta_c : cooling_effective_delta_c(R, p.theta_rec, p.D, p.delta_c, p.thermal)
            s = sin(p.theta_rec)
            c = cos(p.theta_rec)
            geom = g_theta_fast(p.theta_rec) / (s * (2.0 + c))
            du[1] = -(p.D * delta_c_eff / (p.rho_l * R)) * geom
            return nothing
        end

        cb_ca = ContinuousCallback((u, t, integrator) -> u[1] - 1e-9, terminate!)
        prob_ca = ODEProblem(rhs_ca, [cfg.R0], (0.0, t_rem), cfg)
        sol_ca = solve(prob_ca, solver; callback=cb_ca, reltol=rtol, abstol=atol)

        # Concatenate solutions
        t_vec = vcat(sol_cr.t, [t1 + tau for tau in sol_ca.t[2:end]])
        R_vec = vcat(fill(cfg.R0, length(sol_cr.t)), [u[1] for u in sol_ca.u[2:end]])
        theta_vec = vcat([u[1] for u in sol_cr.u], fill(cfg.theta_rec, length(sol_ca.t) - 1))
        V_vec = [sessile_volume_from_r(R_vec[i], theta_vec[i]) for i in eachindex(t_vec)]

        return SessileSolution(t_vec, R_vec, theta_vec, V_vec, :SS, t_ana)

    else # :SJ mode
        t_ana = lifetime_sj(cfg)
        # Numerical realization of periodic pinning/depinning jumps
        t_vec = [0.0]
        R_vec = [cfg.R0]
        theta_vec = [cfg.theta_max]
        V_vec = [sessile_volume_from_r(cfg.R0, cfg.theta_max)]

        R_n = cfg.R0
        t_curr = 0.0
        s_min, c_min = sin(cfg.theta_min), cos(cfg.theta_min)
        s_max, c_max = sin(cfg.theta_max), cos(cfg.theta_max)
        lambda = (s_min / s_max) * ((2.0 + c_min) / (2.0 + c_max)) * (((1.0 + c_max) / (1.0 + c_min))^2)

        for _ in 1:20
            if R_n <= 1e-6 || t_curr >= cfg.t_end
                break
            end
            # Pinning step with evaporative cooling if thermal parameters are provided
            delta_c_eff = cfg.thermal === nothing ? cfg.delta_c : cooling_effective_delta_c(R_n, (cfg.theta_min + cfg.theta_max) / 2.0, cfg.D, cfg.delta_c, cfg.thermal)
            dt_step = (cfg.rho_l * (R_n^2) / (cfg.D * delta_c_eff)) * quadgk(th -> 1.0 / g_theta_fast(th), cfg.theta_min, cfg.theta_max)[1]
            t_curr += dt_step
            push!(t_vec, t_curr)
            push!(R_vec, R_n)
            push!(theta_vec, cfg.theta_min)
            push!(V_vec, sessile_volume_from_r(R_n, cfg.theta_min))

            # Jump step: instantaneous contact line depinning at constant volume
            R_n = R_n * (lambda^(1.0 / 3.0))
            push!(t_vec, t_curr)
            push!(R_vec, R_n)
            push!(theta_vec, cfg.theta_max)
            push!(V_vec, sessile_volume_from_r(R_n, cfg.theta_max))
        end

        return SessileSolution(t_vec, R_vec, theta_vec, V_vec, :SJ, t_ana)
    end
end
