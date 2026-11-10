# Surface Temperature Surrogates: Strategy D (Chebyshev Polynomials & Bicubic LUT)

using LinearAlgebra
using Base.Threads

export ChebyshevTsSurrogate, BicubicLutTsSurrogate
export build_chebyshev_ts_surrogate, eval_chebyshev_ts
export build_lut_ts_surrogate, eval_lut_ts
export get_or_build_ts_surrogate, clear_ts_surrogate_cache!

# -----------------------------------------------------------------------------
# GLOBAL SURROGATE CACHE
# -----------------------------------------------------------------------------
const TS_SURROGATE_CACHE = Dict{Tuple, Any}()
const TS_CACHE_LOCK = ReentrantLock()

"""
    clear_ts_surrogate_cache!()

Clears the global thread-safe cache of precomputed Ts surrogates.
"""
function clear_ts_surrogate_cache!()
    lock(TS_CACHE_LOCK) do
        empty!(TS_SURROGATE_CACHE)
    end
    return nothing
end

# -----------------------------------------------------------------------------
# HELPER: Clenshaw Recurrence for 1D Chebyshev Sum
# -----------------------------------------------------------------------------
@inline function clenshaw_1d(c::AbstractVector{Float64}, x::Float64)
    N = length(c)
    if N == 1
        return c[1]
    end
    b2 = 0.0
    b1 = 0.0
    two_x = 2.0 * x
    @inbounds for k in N:-1:2
        b0 = two_x * b1 - b2 + c[k]
        b2 = b1
        b1 = b0
    end
    @inbounds return x * b1 - b2 + c[1]
end

# -----------------------------------------------------------------------------
# STRATEGY D1: CHEBYSHEV POLYNOMIAL SURROGATE
# -----------------------------------------------------------------------------

"""
    build_chebyshev_ts_surrogate(params::DropletParams; d_max=1500e-6, T_min=nothing, T_max=nothing, n_xi=8, n_T=12)

Constructs a 2D (or 1D when Re=0) Chebyshev polynomial surrogate for Ts(rs, T_bar)
using transformed coordinate xi = sqrt(rs).
"""
function build_chebyshev_ts_surrogate(
    params::DropletParams;
    d_max::Real=1500e-6,
    T_min::Union{Nothing, Real}=nothing,
    T_max::Union{Nothing, Real}=nothing,
    n_xi::Int=8,
    n_T::Int=12
)
    rs_max = Float64(d_max) / 2.0
    xi_max = sqrt(rs_max)

    # Ambient wet bulb estimate for lower temperature bound
    T_inf_val = Float64(params.T_inf)
    twb_est = stull_twb(T_inf_val, max(1.0, Float64(params.Y_inf) * 1000.0))
    t_min = T_min === nothing ? min(twb_est - 5.0, 0.0) : Float64(T_min)
    t_max = T_max === nothing ? max(T_inf_val + 5.0, 50.0) : Float64(T_max)

    is_1d = (!params.terminal_velocity) && (params.U_rel == 0.0)
    N_xi = is_1d ? 1 : max(3, n_xi)
    N_T  = max(4, n_T)

    # Chebyshev-Gauss-Lobatto nodes on [-1, 1]
    # x_k = -cos(pi * k / (N-1)), so k=0 -> -1, k=N-1 -> +1
    nodes_T = [-cos(π * j / (N_T - 1)) for j in 0:(N_T - 1)]
    nodes_xi = is_1d ? [0.0] : [-cos(π * i / (N_xi - 1)) for i in 0:(N_xi - 1)]

    # Physical coordinates
    T_vals = [t_min + 0.5 * (node + 1.0) * (t_max - t_min) for node in nodes_T]
    xi_vals = is_1d ? [xi_max * 0.5] : [0.5 * (node + 1.0) * xi_max for node in nodes_xi]

    # Evaluate exact roots at nodes using warm-started Newton
    F = zeros(Float64, N_xi, N_T)
    # Use Newton or Halley for node solve
    p_solve = DropletParams(
        T_inf = params.T_inf,
        U_rel = params.U_rel,
        Y_inf = params.Y_inf,
        P_atm = params.P_atm,
        profile = params.profile,
        correlation = params.correlation,
        as1989_blowing = params.as1989_blowing,
        fixed_chi_m = params.fixed_chi_m,
        model_variant = params.model_variant,
        terminal_velocity = params.terminal_velocity,
        isothermal = params.isothermal,
        T_d_fixed = params.T_d_fixed,
        thermal_dilation = params.thermal_dilation,
        variable_cp_factor = params.variable_cp_factor,
        spheroid_epsilon = params.spheroid_epsilon,
        zhou_moving_boundary = params.zhou_moving_boundary,
        reference_state = params.reference_state,
        use_wilke_mixtures = params.use_wilke_mixtures,
        ts_solver = :newton
    )

    last_root = (t_min + t_max) * 0.5
    for j in 1:N_T
        tb = T_vals[j]
        rho_l = liquid_density_water(tb + 273.15)
        for i in 1:N_xi
            rs = max(xi_vals[i]^2, 1e-9)
            m = rho_l * (4.0 * π / 3.0) * (rs^3)
            liq = precompute_liquid_state(m, tb, p_solve)
            root = solve_surface_temperature_newton(liq, p_solve; tol=1e-8, maxiter=30)
            F[i, j] = root
            last_root = root
        end
    end

    # DCT Type-I matrix projection to obtain Chebyshev expansion coefficients
    # c_{k, l} = (2 / (N_xi-1)) * (2 / (N_T-1)) * (1 / (d_k * d_l)) * sum( ... )
    function dct1_matrix(N::Int)
        if N == 1
            return ones(Float64, 1, 1)
        end
        M = zeros(Float64, N, N)
        for k in 0:(N - 1)
            d_k = (k == 0 || k == N - 1) ? 2.0 : 1.0
            scale = 2.0 / ((N - 1) * d_k)
            for i in 0:(N - 1)
                d_i = (i == 0 || i == N - 1) ? 2.0 : 1.0
                # Note: node x_i = -cos(pi * i / (N - 1)) = cos(pi * (N - 1 - i) / (N - 1))
                # angle = pi * (N - 1 - i) / (N - 1)
                theta_i = π * (N - 1 - i) / (N - 1)
                M[k + 1, i + 1] = (scale / d_i) * cos(k * theta_i)
            end
        end
        return M
    end

    M_T = dct1_matrix(N_T)
    M_xi = is_1d ? ones(Float64, 1, 1) : dct1_matrix(N_xi)

    # C = M_xi * F * M_T'
    coeffs = M_xi * F * transpose(M_T)

    return ChebyshevTsSurrogate(xi_max, t_min, t_max, is_1d, coeffs)
end

"""
    eval_chebyshev_ts(surr::ChebyshevTsSurrogate, r_s::Real, T_bar::Real) -> Real
"""
function eval_chebyshev_ts(surr::ChebyshevTsSurrogate, r_s::Real, T_bar::Real)
    # ForwardDiff.Dual propagation branch
    if r_s isa ForwardDiff.Dual || T_bar isa ForwardDiff.Dual
        r_val = Float64(ForwardDiff.value(r_s))
        T_val = Float64(ForwardDiff.value(T_bar))
        val = eval_chebyshev_ts(surr, r_val, T_val)

        eps_r = max(1e-12, r_val * 1e-4)
        eps_T = 1e-4
        dTs_dr = (eval_chebyshev_ts(surr, r_val + eps_r, T_val) - eval_chebyshev_ts(surr, max(1e-12, r_val - eps_r), T_val)) / (2.0 * eps_r)
        dTs_dT = (eval_chebyshev_ts(surr, r_val, T_val + eps_T) - eval_chebyshev_ts(surr, r_val, T_val - eps_T)) / (2.0 * eps_T)

        part_r = (r_s isa ForwardDiff.Dual) ? ForwardDiff.partials(r_s) : zero(ForwardDiff.partials(T_bar))
        part_T = (T_bar isa ForwardDiff.Dual) ? ForwardDiff.partials(T_bar) : zero(ForwardDiff.partials(r_s))
        partials = dTs_dr * part_r + dTs_dT * part_T
        return ForwardDiff.Dual{ForwardDiff.tagtype(typeof(r_s isa ForwardDiff.Dual ? r_s : T_bar))}(val, partials)
    end

    # Float64 evaluation
    tb = Float64(T_bar)
    y = 2.0 * (tb - surr.T_min) / (surr.T_max - surr.T_min) - 1.0
    y_c = clamp(y, -1.0, 1.0)

    if surr.is_1d
        c_1d = @view surr.coeffs[1, :]
        return clenshaw_1d(c_1d, y_c)
    end

    rs = max(Float64(r_s), 0.0)
    xi = sqrt(rs)
    x = 2.0 * (xi / surr.xi_max) - 1.0
    x_c = clamp(x, -1.0, 1.0)

    # 2D Clenshaw
    N_xi, N_T = size(surr.coeffs)
    a = zeros(Float64, N_xi)
    @inbounds for i in 1:N_xi
        row = @view surr.coeffs[i, :]
        a[i] = clenshaw_1d(row, y_c)
    end
    return clenshaw_1d(a, x_c)
end

# -----------------------------------------------------------------------------
# STRATEGY D2: PRECOMPUTED 2D BICUBIC HERMITE / CATMULL-ROM LOOKUP TABLE
# -----------------------------------------------------------------------------

"""
    build_lut_ts_surrogate(params::DropletParams; d_max=1500e-6, T_min=nothing, T_max=nothing, n_xi=32, n_T=32)

Constructs an equispaced 2D grid in (xi = sqrt(rs), T_bar) with bicubic Hermite (Catmull-Rom) interpolation.
"""
function build_lut_ts_surrogate(
    params::DropletParams;
    d_max::Real=1500e-6,
    T_min::Union{Nothing, Real}=nothing,
    T_max::Union{Nothing, Real}=nothing,
    n_xi::Int=32,
    n_T::Int=32
)
    rs_max = Float64(d_max) / 2.0
    xi_max = sqrt(rs_max)

    T_inf_val = Float64(params.T_inf)
    twb_est = stull_twb(T_inf_val, max(1.0, Float64(params.Y_inf) * 1000.0))
    t_min = T_min === nothing ? min(twb_est - 5.0, 0.0) : Float64(T_min)
    t_max = T_max === nothing ? max(T_inf_val + 5.0, 50.0) : Float64(T_max)

    N_xi = max(4, n_xi)
    N_T  = max(4, n_T)

    inv_dxi = (N_xi - 1) / xi_max
    inv_dT  = (N_T - 1) / (t_max - t_min)

    xi_vals = range(0.0, xi_max, length=N_xi)
    T_vals  = range(t_min, t_max, length=N_T)

    data = zeros(Float64, N_xi, N_T)

    p_solve = DropletParams(
        T_inf = params.T_inf,
        U_rel = params.U_rel,
        Y_inf = params.Y_inf,
        P_atm = params.P_atm,
        profile = params.profile,
        correlation = params.correlation,
        as1989_blowing = params.as1989_blowing,
        fixed_chi_m = params.fixed_chi_m,
        model_variant = params.model_variant,
        terminal_velocity = params.terminal_velocity,
        isothermal = params.isothermal,
        T_d_fixed = params.T_d_fixed,
        thermal_dilation = params.thermal_dilation,
        variable_cp_factor = params.variable_cp_factor,
        spheroid_epsilon = params.spheroid_epsilon,
        zhou_moving_boundary = params.zhou_moving_boundary,
        reference_state = params.reference_state,
        use_wilke_mixtures = params.use_wilke_mixtures,
        ts_solver = :newton
    )

    for j in 1:N_T
        tb = T_vals[j]
        rho_l = liquid_density_water(tb + 273.15)
        for i in 1:N_xi
            rs = max(xi_vals[i]^2, 1e-9)
            m = rho_l * (4.0 * π / 3.0) * (rs^3)
            liq = precompute_liquid_state(m, tb, p_solve)
            data[i, j] = solve_surface_temperature_newton(liq, p_solve; tol=1e-8, maxiter=30)
        end
    end

    return BicubicLutTsSurrogate(xi_max, t_min, t_max, inv_dxi, inv_dT, data)
end

@inline function catmull_rom(p0::Float64, p1::Float64, p2::Float64, p3::Float64, t::Float64)
    t2 = t * t
    t3 = t2 * t
    return 0.5 * (
        (2.0 * p1) +
        (-p0 + p2) * t +
        (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2 +
        (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3
    )
end

"""
    eval_lut_ts(surr::BicubicLutTsSurrogate, r_s::Real, T_bar::Real) -> Real
"""
function eval_lut_ts(surr::BicubicLutTsSurrogate, r_s::Real, T_bar::Real)
    # ForwardDiff.Dual propagation branch
    if r_s isa ForwardDiff.Dual || T_bar isa ForwardDiff.Dual
        r_val = Float64(ForwardDiff.value(r_s))
        T_val = Float64(ForwardDiff.value(T_bar))
        val = eval_lut_ts(surr, r_val, T_val)

        eps_r = max(1e-12, r_val * 1e-4)
        eps_T = 1e-4
        dTs_dr = (eval_lut_ts(surr, r_val + eps_r, T_val) - eval_lut_ts(surr, max(1e-12, r_val - eps_r), T_val)) / (2.0 * eps_r)
        dTs_dT = (eval_lut_ts(surr, r_val, T_val + eps_T) - eval_lut_ts(surr, r_val, T_val - eps_T)) / (2.0 * eps_T)

        part_r = (r_s isa ForwardDiff.Dual) ? ForwardDiff.partials(r_s) : zero(ForwardDiff.partials(T_bar))
        part_T = (T_bar isa ForwardDiff.Dual) ? ForwardDiff.partials(T_bar) : zero(ForwardDiff.partials(r_s))
        partials = dTs_dr * part_r + dTs_dT * part_T
        return ForwardDiff.Dual{ForwardDiff.tagtype(typeof(r_s isa ForwardDiff.Dual ? r_s : T_bar))}(val, partials)
    end

    N_xi, N_T = size(surr.data)
    rs = max(Float64(r_s), 0.0)
    xi = sqrt(rs)
    tb = Float64(T_bar)

    # Dimension 1: xi
    u = clamp(xi * surr.inv_dxi, 0.0, Float64(N_xi - 1))
    i = clamp(floor(Int, u) + 1, 1, N_xi - 1)
    t_xi = u - (i - 1)

    i0 = max(1, i - 1)
    i1 = i
    i2 = min(N_xi, i + 1)
    i3 = min(N_xi, i + 2)

    # Dimension 2: T_bar
    v = clamp((tb - surr.T_min) * surr.inv_dT, 0.0, Float64(N_T - 1))
    j = clamp(floor(Int, v) + 1, 1, N_T - 1)
    t_T = v - (j - 1)

    j0 = max(1, j - 1)
    j1 = j
    j2 = min(N_T, j + 1)
    j3 = min(N_T, j + 2)

    # Bicubic interpolation
    data = surr.data
    c0 = catmull_rom(data[i0, j0], data[i0, j1], data[i0, j2], data[i0, j3], t_T)
    c1 = catmull_rom(data[i1, j0], data[i1, j1], data[i1, j2], data[i1, j3], t_T)
    c2 = catmull_rom(data[i2, j0], data[i2, j1], data[i2, j2], data[i2, j3], t_T)
    c3 = catmull_rom(data[i3, j0], data[i3, j1], data[i3, j2], data[i3, j3], t_T)

    return catmull_rom(c0, c1, c2, c3, t_xi)
end

# -----------------------------------------------------------------------------
# AUTOMATIC SURROGATE DISPATCH & CACHING
# -----------------------------------------------------------------------------

"""
    get_or_build_ts_surrogate(params::DropletParams, method::Symbol; d_max=1500e-6)

Retrieves a cached surrogate matching `params` and `method` (:chebyshev or :lut),
or builds, caches, and returns a new one if not found.
"""
function get_or_build_ts_surrogate(params::DropletParams, method::Symbol; d_max::Real=1500e-6)
    if params.ts_surrogate !== nothing
        return params.ts_surrogate
    end

    key = (
        Float64(params.T_inf),
        Float64(params.Y_inf),
        Float64(params.P_atm),
        Float64(params.U_rel),
        params.profile,
        params.correlation,
        params.as1989_blowing,
        params.terminal_velocity,
        params.use_wilke_mixtures,
        params.reference_state,
        params.model_variant,
        method,
        Float64(d_max)
    )

    lock(TS_CACHE_LOCK) do
        if haskey(TS_SURROGATE_CACHE, key)
            return TS_SURROGATE_CACHE[key]
        end
        surr = if method === :lut
            build_lut_ts_surrogate(params; d_max=d_max)
        else
            build_chebyshev_ts_surrogate(params; d_max=d_max)
        end
        TS_SURROGATE_CACHE[key] = surr
        return surr
    end
end
