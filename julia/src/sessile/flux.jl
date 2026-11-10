"""
    flux.jl

Global and local evaporative mass fluxes for sessile droplets with substrate cooling.
Based on Popov (2005) and Wilson & D'Ambrosio (2023).
"""

"""
    SessileThermalParams

Substrate and latent heat parameters for conductive cooling correction.
"""
Base.@kwdef struct SessileThermalParams{T<:Real}
    L::T              # Latent heat of vaporization [J/kg]
    h_s::T            # Substrate thickness [m]
    k_s_eff::T        # Effective substrate thermal conductivity [W/(m*K)]
    dcsat_dT::T       # Derivative of saturation vapor concentration w.r.t. temperature [kg/(m^3*K)]
end

"""
    cooling_effective_delta_c(R, theta, D, delta_c0, thermal::SessileThermalParams) -> Real

Effective concentration difference Δc_eff accounting for conductive substrate cooling.
Wilson & D'Ambrosio (2023) / Popov (2005).
"""
function cooling_effective_delta_c(
    R::Real, theta::Real, D::Real, delta_c0::T, thermal::SessileThermalParams
) where {T<:Real}
    th = Float64(theta)
    g_val = g_theta_fast(th)
    c = cos(th)
    f_th = g_val / ((1.0 + c)^2)

    # Dimensionless thermal parameter A = L * D * h_s * (d c_sat / dT) / k_s_eff
    A = thermal.L * D * thermal.h_s * thermal.dcsat_dT / thermal.k_s_eff
    return delta_c0 / (one(T) + T(A * f_th / R))
end

"""
    popov_global_flux(R, theta, D, delta_c; thermal=nothing) -> Real

Total evaporative mass loss rate F [kg/s] from a spherical cap:
    F = π * D * Δc_eff * R * g(θ) / (1 + cos θ)^2
"""
function popov_global_flux(
    R::T, theta::Real, D::Real, delta_c::Real;
    thermal::Union{Nothing, SessileThermalParams}=nothing
) where {T<:Real}
    th = Float64(theta)
    delta_c_eff = thermal === nothing ? T(delta_c) : cooling_effective_delta_c(R, theta, D, T(delta_c), thermal)
    g_val = T(g_theta_fast(th))
    c = cos(T(theta))
    f_th = g_val / ((one(T) + c)^2)
    return T(π) * T(D) * delta_c_eff * R * f_th
end

"""
    estimate_interface_temperature_drop(R, theta, D, delta_c_eff, thermal::SessileThermalParams) -> Real

Mean interface temperature reduction ΔT [K] induced by latent heat of evaporation.
"""
function estimate_interface_temperature_drop(
    R::Real, theta::Real, D::Real, delta_c_eff::Real, thermal::SessileThermalParams
)
    th = Float64(theta)
    g_val = g_theta_fast(th)
    c = cos(th)
    f_th = g_val / ((1.0 + c)^2)
    F = π * R * D * delta_c_eff * f_th
    flux_density = thermal.L * F / (π * (R^2))
    return flux_density * thermal.h_s / thermal.k_s_eff
end

"""
    local_flux_popov_series(r, R, theta, D, delta_c; N=50) -> Real

Local evaporative mass flux J(r) [kg/(m^2*s)] across substrate contact line using Popov (2005) series.
"""
function local_flux_popov_series(r::Real, R::Real, theta::Real, D::Real, delta_c::Real; N::Int=50)
    if r >= R
        return 0.0
    end
    c = cos(theta)
    base = sin(theta) / (1.0 + c)
    s = 0.0
    r_ratio = r / R

    for n in 1:N
        term = (1.0 / (4.0 * n^2 - 1.0)) * (sin(2.0 * n * theta) / ((1.0 + c)^(2.0 * n))) * (r_ratio^(2.0 * n))
        s += term
    end

    return (D * delta_c / R) * (base + (4.0 / π) * s)
end

"""
    local_flux_deegan(r, R, theta, D, delta_c) -> Real

Deegan et al. (1997) near-contact-line singular flux scaling:
    J(r) ∝ (1 - r/R)^(λ - 1), where λ = π / (2(π - θ)).
"""
function local_flux_deegan(r::Real, R::Real, theta::Real, D::Real, delta_c::Real)
    if r >= R
        return 0.0
    end
    th = Float64(theta)
    lambda = th < π ? π / (2.0 * (π - th)) : 0.5
    eta = max(1e-14, 1.0 - (r / R))
    edge = eta^(lambda - 1.0)
    ampl = sin(th) / (1.0 + cos(th))
    return (D * delta_c / R) * edge * ampl
end

# =============================================================================
# PICKNETT & BEXON (1977) CAPACITANCE & GLOBAL FLUX
# =============================================================================

"""
    picknett_bexon_capacitance(theta::Real) -> Real

Calculates the dimensionless capacitance ratio f(θ) = C / R for a spherical cap lens
using the two-branch polynomials from Picknett & Bexon (1977, Eqs. 5a-5b; Wilson & D'Ambrosio 2023 Eq. 11):
    f(θ) = (0.6366*θ + 0.09591*θ^2 - 0.06144*θ^3) / sin(θ)                           for 0 < θ < 0.175 rad (10°)
    f(θ) = (0.00008957 + 0.6333*θ + 0.1160*θ^2 - 0.08878*θ^3 + 0.01033*θ^4) / sin(θ) for 0.175 ≤ θ ≤ π
"""
function picknett_bexon_capacitance(theta::Real)
    th = clamp(Float64(theta), 0.0, π)
    if th < 1e-4
        return 0.6366197723675814 # 2/π limit
    elseif th < 0.175
        s = sin(th)
        num = 0.6366 * th + 0.09591 * (th^2) - 0.06144 * (th^3)
        return num / s
    else
        s = sin(th)
        if s < 1e-4
            return 4.0 * log(2.0) / π
        end
        num = 0.00008957 + 0.6333 * th + 0.1160 * (th^2) - 0.08878 * (th^3) + 0.01033 * (th^4)
        return num / s
    end
end

"""
    picknett_bexon_flux(R::Real, theta::Real, D::Real, delta_c::Real) -> Real

Total evaporative mass flux F [kg/s] from Picknett & Bexon (1977):
    F = 2π * R * D * Δc * f(θ)
"""
function picknett_bexon_flux(R::T, theta::Real, D::Real, delta_c::Real) where {T<:Real}
    f_th = T(picknett_bexon_capacitance(theta))
    return T(2.0 * π) * R * T(D) * T(delta_c) * f_th
end

# =============================================================================
# BOURGES-MONNIER & SHANAHAN (1995) FLUX
# =============================================================================

"""
    bourges_monnier_flux(R::Real, theta::Real, D::Real, delta_c::Real) -> Real

Closed-form approximation for total evaporation flux from Bourges-Monnier & Shanahan (1995)
(Wilson & D'Ambrosio 2023, p. 5):
    F = 2π * D * Δc * R * cot(θ) / (-ln(1 - cos(θ)))
"""
function bourges_monnier_flux(R::T, theta::Real, D::Real, delta_c::Real) where {T<:Real}
    th = clamp(Float64(theta), 1e-4, π - 1e-4)
    c = cos(th)
    cot_th = c / sin(th)
    denom = -log(max(1.0 - c, 1e-12))
    return T(2.0 * π) * T(D) * T(delta_c) * R * T(cot_th / denom)
end

# =============================================================================
# NON-EQUILIBRIUM (NEOS) KINETIC FLUX & GAS CONDUCTION INDEX
# =============================================================================

"""
    local_flux_neos(r::Real, R::Real, theta::Real, K_evap::Real, K_kin::Real) -> Real

One-sided non-equilibrium (NEOS) kinetic local evaporative flux profile
(Burelbach et al. 1988, Murisic & Kondic 2011, Wilson & D'Ambrosio 2023 Eq. 31):
    J(r) = K_evap / (h(r) + K_kin)
where h(r) = (θ*R/2) * (1 - (r/R)^2).
"""
function local_flux_neos(r::T, R::Real, theta::Real, K_evap::Real, K_kin::Real) where {T<:Real}
    r_val = clamp(r, zero(T), T(R))
    r_ratio = r_val / T(R)
    h_local = (T(theta) * T(R) / T(2.0)) * (one(T) - r_ratio^2)
    return T(K_evap) / (h_local + T(K_kin))
end

"""
    gas_conduction_index(T_inf::Real, T_s::Real, J_mean::Real, L_v::Real, k_g::Real, R::Real) -> Real

Dimensionless gas heat conduction index HC_g from Xu & Xu (2025):
    HC_g = k_g * (T_inf - T_s) / (J_mean * L_v * R)
"""
function gas_conduction_index(
    T_inf::T, T_s::Real, J_mean::Real, L_v::Real, k_g::Real, R::Real
) where {T<:Real}
    delta_T = max(T_inf - T(T_s), zero(T))
    denom = max(T(J_mean) * T(L_v) * T(R), T(1e-12))
    return T(k_g) * delta_T / denom
end

# =============================================================================
# WRAY ET AL. (2020) COOPERATIVE MULTI-DROPLET VAPOR SHIELDING SYSTEM
# =============================================================================

"""
    cooperative_shielding_fluxes(radii::Vector{<:Real}, positions::Matrix{<:Real}, D::Real, delta_c::Real, thetas::Vector{<:Real}) -> Vector{Float64}

Solves the N x N linear algebraic system for cooperative vapor shielding fluxes among
an arbitrary array of N sessile droplets (Wray et al. 2020, Masoud et al. 2021, Wilson & D'Ambrosio 2023):
    F_m + (2/π) * ∑_{k ≠ m} F_k * arcsin(R_k / d_mk) = F_{m,iso}
where `positions` is an N x 2 matrix of (x, y) coordinates and `F_{m,iso}` is the isolated Popov flux.
"""
function cooperative_shielding_fluxes(
    radii::Vector{<:Real}, positions::AbstractMatrix{<:Real},
    D::Real, delta_c::Real, thetas::Vector{<:Real}
)
    N = length(radii)
    @assert size(positions, 1) == N && size(positions, 2) >= 2 "positions must be N x 2"
    @assert length(thetas) == N "thetas must have length N"

    M = Matrix{Float64}(I, N, N)
    F_iso = zeros(Float64, N)

    for i in 1:N
        R_i = Float64(radii[i])
        th_i = Float64(thetas[i])
        F_iso[i] = popov_global_flux(R_i, th_i, D, delta_c)

        for j in 1:N
            if i != j
                dx = Float64(positions[i, 1] - positions[j, 1])
                dy = Float64(positions[i, 2] - positions[j, 2])
                dist = sqrt(dx^2 + dy^2)
                R_j = Float64(radii[j])
                arg = clamp(R_j / max(dist, R_j + 1e-12), 0.0, 1.0)
                M[i, j] = (2.0 / π) * asin(arg)
            end
        end
    end

    # Solve linear system M * F = F_iso
    return M \ F_iso
end
