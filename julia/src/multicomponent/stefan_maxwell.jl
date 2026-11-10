"""
    Stefan-Maxwell exact analytical matrix exponential and Stefan-Fuchs radial convection models.
    Based on Tonini & Cossali (2015) and Ordoubadi et al. (2019).
"""

using LinearAlgebra

"""
    stefan_maxwell_matrices(N_tot_hat, mu, omega)

Construct the Tonini & Cossali (2015) Stefan-Maxwell coefficient matrix `A` and source vector `B`.
- `N_tot_hat`: Dimensionless total molar evaporation rate
- `mu`: Vector of species molar evaporation fractions (length n)
- `omega`: Matrix of relative binary diffusivities omega[p, k] = D_ref / D[p, k] (size n x (n+1), where column 1:n are species and column n+1 or index 0 is carrier gas)
"""
function stefan_maxwell_matrices(N_tot_hat::Real, mu::AbstractVector{<:Real}, omega_species::AbstractMatrix{<:Real}, omega_carrier::AbstractVector{<:Real})
    n = length(mu)
    A = zeros(Float64, n, n)
    B = zeros(Float64, n)

    for p in 1:n
        B[p] = N_tot_hat * mu[p] * omega_carrier[p]
        for j in 1:n
            if p == j
                sum_k = sum(omega_species[p, k] * mu[k] for k in 1:n if k != p; init=0.0) + omega_carrier[p] * mu[p]
                A[p, p] = -N_tot_hat * sum_k
            else
                A[p, j] = N_tot_hat * mu[p] * (omega_species[p, j] - omega_carrier[p])
            end
        end
    end
    return A, B
end

"""
    stefan_maxwell_residual(N_tot_hat, mu, W_s, W_inf, omega_species, omega_carrier)

Evaluate the analytical matrix exponential residual vector:
```math
\\mathbf{R} = e^{\\mathbf{A}} (\\mathbf{A} \\mathbf{W}_\\infty + \\mathbf{B}) - (\\mathbf{A} \\mathbf{W}_s + \\mathbf{B}) = \\mathbf{0}
```
"""
function stefan_maxwell_residual(N_tot_hat::Real, mu::AbstractVector{<:Real},
                                 W_s::AbstractVector{<:Real}, W_inf::AbstractVector{<:Real},
                                 omega_species::AbstractMatrix{<:Real}, omega_carrier::AbstractVector{<:Real})
    A, B = stefan_maxwell_matrices(N_tot_hat, mu, omega_species, omega_carrier)
    expA = exp(A)
    res = expA * (A * W_inf + B) - (A * W_s + B)
    return res
end

"""
    stefan_fuchs_multicomponent_fluxes(d, rho_t, D_inert_m, y_s_inert, y_inf_inert, y_s, y_inf, D_species)

Calculate competitive individual solvent mass evaporation rates [kg/s] according to Ordoubadi et al. (2019):
```math
\\dot{m}_t = -2 \\pi d \\rho_t \\mathcal{D}_{\\text{inert},m} \\ln\\left( \\frac{y_{\\infty,\\text{inert}}}{y_{s,\\text{inert}}} \\right)
```
```math
\\dot{m}_i = \\dot{m}_t \\left[ \\frac{y_{\\infty,i} \\beta_i - y_{s,i}}{\\beta_i - 1} \\right], \\qquad \\beta_i = \\left( \\frac{y_{s,\\text{inert}}}{y_{\\infty,\\text{inert}}} \\right)^{\\frac{\\mathcal{D}_{\\text{inert},m}}{\\mathcal{D}_{i,m}}}
```
"""
function stefan_fuchs_multicomponent_fluxes(d::Real, rho_t::Real, D_inert_m::Real,
                                            y_s_inert::Real, y_inf_inert::Real,
                                            y_s::AbstractVector{<:Real}, y_inf::AbstractVector{<:Real},
                                            D_species::AbstractVector{<:Real})
    if y_s_inert <= 0.0 || y_inf_inert <= 0.0 || d <= 0.0
        return 0.0, zeros(Float64, length(y_s))
    end

    # Total mass evaporation rate
    m_dot_t = 2.0 * pi * float(d) * float(rho_t) * float(D_inert_m) * log(float(y_inf_inert) / float(y_s_inert))

    # Individual rates
    m_dot_i = zeros(Float64, length(y_s))
    inert_ratio = float(y_s_inert) / float(y_inf_inert)

    for i in 1:length(y_s)
        diff_exponent = float(D_inert_m) / float(D_species[i])
        beta_i = inert_ratio^diff_exponent

        if abs(beta_i - 1.0) < 1e-7
            # Exact L'Hôpital limit when beta_i -> 1 (diffusivities equal or dilute limit m_dot_t -> 0)
            m_dot_i[i] = 2.0 * π * float(d) * float(rho_t) * float(D_species[i]) * (y_s[i] - y_inf[i])
        else
            m_dot_i[i] = m_dot_t * (y_inf[i] * beta_i - y_s[i]) / (beta_i - 1.0)
        end
        # Enforce physical positivity for evaporation into low ambient concentrations
        m_dot_i[i] = max(0.0, m_dot_i[i])
    end

    return m_dot_t, m_dot_i
end
