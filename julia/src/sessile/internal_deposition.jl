"""
    internal_deposition.jl

Internal velocity fields, streamlines, and coffee-ring mass deposition
for pinned sessile droplets.
Wilson & D'Ambrosio (2023) Eqs. (24)–(26) and Deegan et al. (1997).
"""

"""
    depth_averaged_radial_velocity(r, R0, theta, D, delta_c, rho_l) -> Real

Depth-averaged radial outward flow u_bar(r) driving the coffee-ring effect.
Wilson & D'Ambrosio (2023) Eq. (24).
"""
function depth_averaged_radial_velocity(
    r::Real, R0::Real, theta::Real, D::Real, delta_c::Real, rho_l::Real
)
    if r <= 1e-15
        return 0.0
    end
    safe_r = min(r, R0 * (1.0 - 1e-12))
    r_ratio_sq = (safe_r / R0)^2

    term1 = 4.0 * D * delta_c / (π * rho_l * theta * safe_r)
    term2 = (1.0 - (1.0 - r_ratio_sq)^1.5) / sqrt(1.0 - r_ratio_sq)
    return term1 * term2
end

"""
    velocity_field_u(r, z, R0, theta, D, delta_c, rho_l) -> Real

Exact radial velocity profile u(r, z) with no-slip at substrate (z=0) and stress-free at free surface.
Wilson & D'Ambrosio (2023) Eq. (25).
"""
function velocity_field_u(
    r::Real, z::Real, R0::Real, theta::Real, D::Real, delta_c::Real, rho_l::Real
)
    if r <= 1e-15 || z <= 0.0
        return 0.0
    end
    safe_r = min(r, R0 * (1.0 - 1e-12))
    r_ratio_sq = (safe_r / R0)^2

    pre = 24.0 * D * delta_c * z / (π * rho_l * (theta^3) * (R0^2) * safe_r)
    h_film = theta * R0 * (1.0 - r_ratio_sq)
    term_bracket = h_film - z
    shape_factor = (1.0 - (1.0 - r_ratio_sq)^1.5) / ((1.0 - r_ratio_sq)^2.5)

    return pre * term_bracket * shape_factor
end

"""
    velocity_field_w(r, z, R0, theta, D, delta_c, rho_l) -> Real

Vertical velocity profile w(r, z) satisfying mass continuity:
    (1/r) ∂(r*u)/∂r + ∂w/∂z = 0
Wilson & D'Ambrosio (2023) Eq. (26).
"""
function velocity_field_w(
    r::Real, z::Real, R0::Real, theta::Real, D::Real, delta_c::Real, rho_l::Real
)
    if z <= 0.0
        return 0.0
    end
    safe_r = min(r, R0 * (1.0 - 1e-12))
    r_ratio_sq = (safe_r / R0)^2

    pre = -4.0 * D * delta_c * (z^2) / (π * rho_l * (theta^3) * (R0^4))
    bracket = 9.0 * theta * R0 * (1.0 - r_ratio_sq) + (4.0 * (1.0 - r_ratio_sq)^1.5 - 10.0) * z
    shape_factor = (1.0 - r_ratio_sq)^(-3.5)

    return pre * bracket * shape_factor
end

"""
    coffee_ring_deposit_scaling(r, R0, theta0) -> Real

Spatial mass deposit profile near the contact line:
    M(r) ∝ (1 - r/R0)^(λ - 1), where λ = π / (2(π - θ₀)).
"""
function coffee_ring_deposit_scaling(r::Real, R0::Real, theta0::Real)
    if r >= R0
        return 0.0
    end
    th = Float64(theta0)
    lambda = th < π ? π / (2.0 * (π - th)) : 0.5
    eta = max(1e-14, 1.0 - (r / R0))
    return eta^(lambda - 1.0)
end

# =============================================================================
# GELDERBLOM ET AL. (2022) UNPINNED (CCA) INTERNAL RADIAL FLOW
# =============================================================================

"""
    depth_averaged_radial_velocity_cca(r, R0, theta, D, delta_c, rho_l) -> Real

Depth-averaged radial flow u_bar_CCA(r) for unpinned Constant Contact Angle (CCA) droplets
(Gelderblom et al. 2022, Eq. 3.2.2). Unlike pinned CR mode, velocity remains finite
and relative velocity matches the receding contact line.
"""
function depth_averaged_radial_velocity_cca(
    r::Real, R0::Real, theta::Real, D::Real, delta_c::Real, rho_l::Real
)
    if r <= 1e-15
        return 0.0
    end
    safe_r = min(r, R0 * (1.0 - 1e-12))
    r_ratio = safe_r / R0
    r_ratio_sq = r_ratio^2

    term1 = 4.0 * D * delta_c / (π * rho_l * theta * safe_r)
    geom = (1.0 - (1.0 - r_ratio_sq)^1.5) / sqrt(1.0 - r_ratio_sq) - 1.5 * r_ratio_sq
    return term1 * geom
end

# =============================================================================
# MARANGONI & THERMAL INVERSION ANALYSIS (Rocha et al. 2024, Wilson 2023)
# =============================================================================

"""
    marangoni_thermal_solutal_analysis(R, DeltaT, Gamma0, mu_l, kappa_l, dsigma_dT, dsigma_dGamma)

Calculates thermal and solutal Marangoni numbers (Ma_T, Ma_Γ) and determines whether
trace surfactant contamination halts or inverts thermal Marangoni flow into buoyancy-dominated
Rayleigh flow (Rocha et al. 2024, JFM):
    Ma_T = R * |∂_T σ| * ΔT / (μ_l * κ_l)
    Ma_Γ = R * |∂_Γ σ| * Γ₀ / (μ_l * κ_l)
Returns `(Ma_T, Ma_Gamma, flow_inverted, ratio)`.
"""
function marangoni_thermal_solutal_analysis(
    R::Real, DeltaT::Real, Gamma0::Real,
    mu_l::Real, kappa_l::Real, dsigma_dT::Real, dsigma_dGamma::Real
)
    denom = max(mu_l * kappa_l, 1e-16)
    Ma_T = R * abs(dsigma_dT) * abs(DeltaT) / denom
    Ma_Gamma = R * abs(dsigma_dGamma) * abs(Gamma0) / denom

    ratio = Ma_T > 0.0 ? Ma_Gamma / Ma_T : Inf
    # Rocha et al. (2024) prove flow inversion occurs at Ma_Γ / Ma_T ≈ 0.8 - 1.0 (0.32% - 0.5% Δσ)
    flow_inverted = ratio >= 0.8

    return (Ma_T = Ma_T, Ma_Gamma = Ma_Gamma, flow_inverted = flow_inverted, ratio = ratio)
end
