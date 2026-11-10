"""
    Evaporacion.Common.Geometry

Spherical and sessile droplet geometric formulations, exact volume derivatives,
and contact angle correction functions g(θ) from Wilson & D'Ambrosio (2023).
"""
module Geometry

using QuadGK

export sphere_surface_from_d, sphere_volume_from_d, sphere_diameter_from_volume,
       sessile_volume_from_r, sessile_surface_from_r, r_from_volume_sessile,
       dV_dtheta, dV_dR,
       g_theta_exact, g_theta_fast, g_theta_hu_larson, g_theta_hu_2014

# =============================================================================
# SPHERICAL DROPLET GEOMETRY
# =============================================================================

sphere_surface_from_d(d::T) where {T<:Real} = T(π) * (d^2)
sphere_volume_from_d(d::T) where {T<:Real} = (T(π) / T(6.0)) * (d^3)
sphere_diameter_from_volume(V::T) where {T<:Real} = cbrt(T(6.0) * max(V, zero(T)) / T(π))

# =============================================================================
# SESSILE DROPLET (SPHERICAL CAP) GEOMETRY
# =============================================================================

"""
    sessile_volume_from_r(R::Real, theta::Real) -> Real

Volume of a spherical cap with contact radius R and contact angle theta (radians).
    V(R, θ) = (π * R^3 / 3) * sin(θ) * (2 + cos(θ)) / (1 + cos(θ))^2
"""
function sessile_volume_from_r(R::T, theta::Real) where {T<:Real}
    th = T(theta)
    s = sin(th)
    c = cos(th)
    factor = s * (T(2.0) + c) / ((one(T) + c)^2)
    return (T(π) * (R^3) / T(3.0)) * factor
end

"""
    sessile_surface_from_r(R::Real, theta::Real) -> Real

Surface area of the liquid-gas interface of a spherical cap.
    A(R, θ) = 2π * R^2 / (1 + cos(θ))
"""
function sessile_surface_from_r(R::T, theta::Real) where {T<:Real}
    th = T(theta)
    c = cos(th)
    return T(2.0 * π) * (R^2) / (one(T) + c)
end

"""
    r_from_volume_sessile(V::Real, theta::Real) -> Real

Contact radius R given volume V and contact angle theta.
"""
function r_from_volume_sessile(V::T, theta::Real) where {T<:Real}
    th = T(theta)
    s = sin(th)
    c = cos(th)
    inv_factor = ((one(T) + c)^2) / (s * (T(2.0) + c))
    return cbrt(T(3.0) * max(V, zero(T)) / T(π) * inv_factor)
end

"""
    dV_dtheta(R::Real, theta::Real) -> Real

Exact partial derivative ∂V/∂θ at constant contact radius R:
    ∂V/∂θ = π * R^3 / (1 + cos(θ))^2
"""
function dV_dtheta(R::T, theta::Real) where {T<:Real}
    th = T(theta)
    c = cos(th)
    den = (one(T) + c)^2
    return T(π) * (R^3) / den
end

"""
    dV_dR(R::Real, theta::Real) -> Real

Exact partial derivative ∂V/∂R at constant contact angle theta:
    ∂V/∂R = π * R^2 * sin(θ) * (2 + cos(θ)) / (1 + cos(θ))^2
"""
function dV_dR(R::T, theta::Real) where {T<:Real}
    th = T(theta)
    s = sin(th)
    c = cos(th)
    return T(π) * (R^2) * s * (T(2.0) + c) / ((one(T) + c)^2)
end

# =============================================================================
# WILSON & D'AMBROSIO (2023) g(θ) FLUX INTEGRAL & LOOKUP
# =============================================================================

"""
    g_theta_exact(theta::Real) -> Float64

Evaluates the exact Wilson & D'Ambrosio (2023) Eq. (10) integral:
    g(θ) = (1 + cos θ)^2 [ tan(θ/2) + 8 ∫_0^∞ cosh^2(θ τ) tanh((π - θ) τ) / sinh(2π τ) dτ ]
"""
function g_theta_exact(theta::Real)
    th = Float64(theta)
    if th <= 1e-6
        return 16.0 / π
    elseif th >= π - 1e-6
        return 0.0
    elseif abs(th - 0.5 * π) < 1e-8
        return 2.0
    end

    integrand = function (tau)
        if tau < 1e-9
            return (π - th) / (2.0 * π)
        end
        # Numerically stable evaluation avoiding overflow in cosh / sinh for large tau:
        # cosh^2(th * tau) / sinh(2π * tau) = (exp(th*tau) + exp(-th*tau))^2 / (2 * sinh(2π*tau))
        # = exp(2*th*tau - 2π*tau) * (1 + exp(-2*th*tau))^2 / (2 * (1 - exp(-4π*tau)))
        a = th * tau
        c = 2.0 * π * tau
        num = (1.0 + exp(-2.0 * a))^2
        den = 2.0 * (1.0 - exp(-2.0 * c))
        kernel = exp(2.0 * a - c) * num / den
        return kernel * tanh((π - th) * tau)
    end

    integral_val, _ = quadgk(integrand, 0.0, 50.0; rtol=1e-8)
    prefactor = (1.0 + cos(th))^2
    return prefactor * (tan(th / 2.0) + 8.0 * integral_val)
end

# Precomputed fast table for g(theta)
const _N_TABLE = 1000
const _THETA_GRID = range(1e-4, π - 1e-4, length=_N_TABLE)
const _G_TABLE = Float64[g_theta_exact(th) for th in _THETA_GRID]

"""
    g_theta_fast(theta::Real) -> Float64

Fast O(1) linear interpolation of g(θ) using precomputed table.
"""
function g_theta_fast(theta::Real)
    th = Float64(theta)
    if th <= 1e-4
        return 16.0 / π
    elseif th >= π - 1e-4
        return 0.0
    end
    # Binary search or index calculation
    dth = (π - 2e-4) / (_N_TABLE - 1)
    idx = clamp(floor(Int, (th - 1e-4) / dth) + 1, 1, _N_TABLE - 1)
    t = (th - _THETA_GRID[idx]) / dth
    return (1.0 - t) * _G_TABLE[idx] + t * _G_TABLE[idx + 1]
end

"""
    g_theta_hu_larson(theta::Real) -> Real

Hu-Larson (2002) empirical approximation:
    g(θ) = (0.27 * θ^2 + 1.30) * (1 + cos θ)^2
"""
function g_theta_hu_larson(theta::T) where {T<:Real}
    c = cos(theta)
    return (T(0.27) * (theta^2) + T(1.30)) * ((one(T) + c)^2)
end

"""
    g_theta_hu_2014(theta::Real) -> Real

Hu et al. (2014) approximation:
    g(θ) = 2 * (1 + cos θ)^1.5
"""
function g_theta_hu_2014(theta::T) where {T<:Real}
    c = cos(theta)
    return T(2.0) * ((one(T) + c)^T(1.5))
end

end # module Geometry
