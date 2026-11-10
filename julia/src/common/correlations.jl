"""
    Evaporacion.Common.Correlations

Dimensionless numbers, transfer correlations (Frössling, Clift), Abramzon-Sirignano (1989)
blowing closures, and Kronig-Brink internal droplet circulation enhancements.
Generic on `T <: Real` for complete ForwardDiff compatibility.
"""
module Correlations

using ..Constants: GRAVITY
using ..Properties: liquid_density_water, gas_density_air, gas_viscosity_air, stull_twb

export calc_re, calc_pr, calc_sc,
       nusselt_froessling, sherwood_froessling, clift_correlation,
       log1p_over_x, blowing_factor, modified_transfer_number,
       TransferSolution, solve_as1989_transfer,
       internal_circulation_enhancement,
       calc_terminal_velocity,
       spheroid_deformation_evaporation_ratio,
       holterman_evaporation_rate, holterman_relaxation_time,
       renaudo_evaporation_parameter, renaudo_extinction_diameter,
       renaudo_deposition_time, renaudo_landing_diameter,
       chiu_group_evaporation_number, chiu_evaporation_regime

# =============================================================================
# DIMENSIONLESS NUMBERS
# =============================================================================

function calc_re(d_m::T, U_ms::Real, rho_g::Real, mu_g::Real) where {T<:Real}
    if mu_g <= 0
        return zero(T)
    end
    return T(rho_g) * abs(T(U_ms)) * d_m / T(mu_g)
end

function calc_pr(mu_g::T, cp_g::Real, k_g::Real) where {T<:Real}
    if k_g <= 0
        return zero(T)
    end
    return mu_g * T(cp_g) / T(k_g)
end

function calc_sc(mu_g::T, rho_g::Real, Dv::Real) where {T<:Real}
    denom = T(rho_g) * T(Dv)
    if denom <= 0
        return zero(T)
    end
    return mu_g / denom
end

# =============================================================================
# BASE NUSSELT & SHERWOOD CORRELATIONS
# =============================================================================

function nusselt_froessling(Re::T, Pr::Real; coef::Real=0.552) where {T<:Real}
    re_pos = max(Re, zero(T))
    pr_pos = max(T(Pr), T(1e-12))
    return T(2.0) + T(coef) * sqrt(re_pos) * (pr_pos^T(1.0 / 3.0))
end

function sherwood_froessling(Re::T, Sc::Real; coef::Real=0.552) where {T<:Real}
    re_pos = max(Re, zero(T))
    sc_pos = max(T(Sc), T(1e-12))
    return T(2.0) + T(coef) * sqrt(re_pos) * (sc_pos^T(1.0 / 3.0))
end

function clift_correlation(Re::T, Pr::Real, Sc::Real) where {T<:Real}
    re_pos = max(Re, zero(T))
    fRe = re_pos <= T(10.0) ? one(T) : (re_pos^T(0.077))
    Nu0 = one(T) + ((one(T) + re_pos * T(Pr))^T(1.0 / 3.0)) * fRe
    Sh0 = one(T) + ((one(T) + re_pos * T(Sc))^T(1.0 / 3.0)) * fRe
    return Nu0, Sh0
end

# =============================================================================
# ABRAMZON-SIRIGNANO (1989) FILM BLOWING FACTOR & CLOSURES
# =============================================================================

"""
    log1p_over_x(x::Real) -> Real

Compute ln(1+x)/x with series expansion near x=0 to avoid loss of precision.
"""
function log1p_over_x(x::T) where {T<:Real}
    if abs(x) < T(1e-7)
        return one(T) - T(0.5) * x + (x^2) / T(3.0)
    end
    return log1p(x) / x
end

"""
    blowing_factor(B::Real) -> Real

Universal Abramzon-Sirignano blowing function:
    F(B) = (1 + B)^0.7 * ln(1 + B) / B
"""
function blowing_factor(B::T) where {T<:Real}
    b = max(B, T(-0.99999999))
    return ((one(T) + b)^T(0.7)) * log1p_over_x(b)
end

"""
    modified_transfer_number(N0::Real, F::Real) -> Real

Modified transfer number with blowing correction:
    N* = 2 + (N0 - 2) / F
"""
function modified_transfer_number(N0::T, F::Real) where {T<:Real}
    f_safe = max(T(F), T(1e-12))
    return T(2.0) + (N0 - T(2.0)) / f_safe
end

struct TransferSolution{T<:Real}
    Nu0::T
    Sh0::T
    Nu_star::T
    Sh_star::T
    B_M::T
    B_T::T
    F_M::T
    F_T::T
    phi::T
end

"""
    solve_as1989_transfer(Re, Pr, Sc, B_M, cp_v, cp_g, Le; coef=0.552, max_iter=100, tol=1e-10)

Solve coupled Abramzon-Sirignano (1989) Eqs. (20, 22, 23) transfer closures iteratively.
"""
function solve_as1989_transfer(
    Re::T, Pr::Real, Sc::Real, B_M::Real,
    cp_v::Real, cp_g::Real, Le::Real;
    coef::Real=0.552,
    max_iter::Int=100,
    tol::Real=1e-10
) where {T<:Real}
    Nu0 = nusselt_froessling(Re, Pr; coef=coef)
    Sh0 = sherwood_froessling(Re, Sc; coef=coef)

    bm = max(T(B_M), T(-0.99999999))
    F_M = blowing_factor(bm)
    Sh_star = modified_transfer_number(Sh0, F_M)

    # Fixed-point iteration for B_T
    bt = bm
    Nu_star = modified_transfer_number(Nu0, one(T))
    F_T = one(T)
    phi = zero(T)

    cp_ratio = T(cp_v) / max(T(cp_g), T(1e-12))
    inv_Le = one(T) / max(T(Le), T(1e-12))

    for _ in 1:max_iter
        F_T = blowing_factor(bt)
        Nu_star = modified_transfer_number(Nu0, F_T)
        phi = cp_ratio * (Sh_star / max(Nu_star, T(1e-12))) * inv_Le
        bt_new = ((one(T) + bm)^phi) - one(T)
        bt_new = max(bt_new, T(-0.99999999))
        if abs(bt_new - bt) < T(tol)
            bt = bt_new
            break
        end
        bt = bt_new
    end

    return TransferSolution(Nu0, Sh0, Nu_star, Sh_star, bm, bt, F_M, F_T, phi)
end

# =============================================================================
# KRONIG-BRINK INTERNAL CIRCULATION ENHANCEMENT
# =============================================================================

"""
    internal_circulation_enhancement(profile, Re, rho_g, mu_l, r_s, U_rel, rho_l, Dab, k_l, cp_l, B_M)

Calculate thermal and mass circulation enhancement factors (chi_l, chi_m, b_q, b_m, U_s).
Uses exact Abramzon & Sirignano (1989) Eq. (30 & 32):
    C_F = 12.69 / (Re^(2/3) * (1 + B_M))
    U_s = (1/32) * U_rel * Re * (mu_g / mu_l) * C_F = rho_g * U_rel^2 * r_s * C_F / (16 * mu_l)
"""
function internal_circulation_enhancement(
    profile::Symbol,
    Re::T, rho_g::Real, mu_l::Real, r_s::Real, U_rel::Real,
    rho_l::Real, Dab::Real, k_l::Real, cp_l::Real, B_M::Real
) where {T<:Real}
    if profile === :uniforme
        return (chi_l = one(T), chi_m = one(T), b_q = zero(T), b_m = zero(T), U_s = zero(T))
    end

    # Abramzon & Sirignano (1989) Eq. (32)
    bm_pos = max(T(B_M), zero(T))
    re_pos = max(Re, zero(T))
    
    C_F = re_pos > T(1e-3) ? T(12.69) / ((re_pos^T(2.0 / 3.0)) * (one(T) + bm_pos)) : zero(T)

    # Surface circulation velocity U_s (Eq. 30)
    U_s = T(rho_g) * (T(U_rel)^2) * T(r_s) * C_F / (T(16.0) * max(T(mu_l), T(1e-12)))
    U_s = max(U_s, T(1e-20))

    # Thermal circulation enhancement chi_l (Abramzon & Sirignano 1989 Eq. 39)
    alpha_l = T(k_l) / (T(rho_l) * T(cp_l))
    Pe_Tl = U_s * T(2.0) * T(r_s) / max(alpha_l, T(1e-16))
    pe_scaled = max(Pe_Tl / T(30.0), T(1e-10))
    chi_l = T(1.86) + T(0.86) * tanh(T(2.245) * log10(pe_scaled))

    # Mass circulation enhancement chi_m (Abramzon & Sirignano 1989 Eq. 39)
    Pe_ml = U_s * T(2.0) * T(r_s) / max(T(Dab), T(1e-16))
    pe_m_scaled = max(Pe_ml / T(30.0), T(1e-10))
    chi_m = T(1.86) + T(0.86) * tanh(T(2.245) * log10(pe_m_scaled))

    if profile === :parabolico
        b_q = T(5.0) * chi_l
        b_m = T(5.0) * chi_m
    elseif profile === :cubico || profile === :zhou_cubic
        # Zhou et al. (2017) Eq. (17): b_q = 1 / 0.086 ≈ 11.6279
        b_q = T(1.0 / 0.086)
        b_m = T(1.0 / 0.086)
    else
        b_q = zero(T)
        b_m = zero(T)
    end

    return (chi_l = chi_l, chi_m = chi_m, b_q = b_q, b_m = b_m, U_s = U_s)
end

"""
    calc_terminal_velocity(d::Real, T_K::Real, P_atm::Real=101325.0; rho_l=nothing, T_liq_K=nothing) -> Real

Calculates the free-fall terminal velocity of a water droplet of diameter d in air,
using the Abramzon & Sirignano (1989) drag correlation C_D = (24/Re)*(1 + Re^(2/3)/6).
Gas properties (rho_g, mu_g) are evaluated at gas film temperature T_K [K].
Liquid droplet density rho_l is evaluated at the droplet bulk temperature T_liq_K [K] (or supplied via `rho_l`).
"""
function calc_terminal_velocity(
    d::T, T_K::Real, P_atm::Real=101325.0;
    rho_l::Union{Nothing, Real}=nothing,
    T_liq_K::Union{Nothing, Real}=nothing
) where {T<:Real}
    if d <= 0
        return zero(T)
    end
    rl = if rho_l !== nothing
        T(rho_l)
    elseif T_liq_K !== nothing
        liquid_density_water(T(T_liq_K))
    else
        liquid_density_water(T(T_K))
    end
    rho_g = gas_density_air(T(T_K), T(P_atm))
    mu_g  = gas_viscosity_air(T(T_K))
    g     = T(GRAVITY)

    v_stokes = g * (d^2) * (rl - rho_g) / (T(18.0) * mu_g)
    v = v_stokes
    for _ in 1:25
        Re = rho_g * v * d / max(mu_g, T(1e-12))
        f = v * (one(T) + (Re^T(2.0 / 3.0)) / T(6.0)) - v_stokes
        df = one(T) + (T(5.0 / 18.0)) * (Re^T(2.0 / 3.0))
        v_next = v - f / df
        if abs(v_next - v) < T(1e-7)
            break
        end
        v = max(v_next, T(1e-8))
    end
    return v
end

"""
    spheroid_deformation_evaporation_ratio(epsilon::Real) -> Real

Calculates the Tonini & Cossali (2015) and Sazhin (2021, Slide 85) ratio of vaporization
rate for an oblate spheroidal droplet (semi-axes ratio a/b = 1 + ε) compared to an
equivalent-volume spherical droplet:
    m_dot_spheroid / m_dot_sphere = sqrt(ε*(2+ε)) / ((1+ε)^(1/3) * atan(sqrt(ε*(2+ε))))
For ε -> 0 (sphere), the ratio smoothly converges to 1 via Taylor series.
"""
function spheroid_deformation_evaporation_ratio(epsilon::T) where {T<:Real}
    eps_val = max(epsilon, zero(T))
    if eps_val < T(1e-4)
        # Taylor expansion around ε = 0: 1 + (2/45)*ε^2 - (2/45)*ε^3 + ...
        return one(T) + (T(2.0 / 45.0)) * (eps_val^2) - (T(2.0 / 45.0)) * (eps_val^3)
    end
    arg = sqrt(eps_val * (T(2.0) + eps_val))
    den = ((one(T) + eps_val)^T(1.0 / 3.0)) * atan(arg)
    return arg / den
end

"""
    holterman_evaporation_rate(D::Real, v_rel::Real, a::Real, b::Real) -> Real

Computes the diameter rate of change dD/dt [m/s] from Holterman (2003, Eq. 3.2):
    dD/dt = -(a / D) * (1 + b * sqrt(D * v_rel))
"""
function holterman_evaporation_rate(D::T, v_rel::Real, a::Real, b::Real) where {T<:Real}
    d_safe = max(D, T(1e-9))
    v_safe = max(T(v_rel), zero(T))
    return -(T(a) / d_safe) * (one(T) + T(b) * sqrt(d_safe * v_safe))
end

"""
    holterman_relaxation_time(D::Real, rho_d::Real, ch::Real, lambda_a::Real, D_vap::Real, M_vap::Real, Lv::Real, T_avg::Real, gamma_psych::Real) -> Real

Thermal relaxation time scale τ_c [s] for a droplet to adjust to wet-bulb equilibrium
from Holterman (2003, Eq. 3.12).
"""
function holterman_relaxation_time(
    D::T, rho_d::Real, ch::Real, lambda_a::Real,
    D_vap::Real, M_vap::Real, Lv::Real, T_avg_K::Real, gamma_psych::Real
) where {T<:Real}
    R_u = T(8314.46) # J/(kmol * K)
    num = T(rho_d) * T(ch) * (D^2)
    den = T(12.0) * (T(lambda_a) + (T(D_vap) * T(M_vap) * T(Lv) / (R_u * T(T_avg_K))) * T(gamma_psych))
    return num / max(den, T(1e-12))
end

# =============================================================================
# RENAUDO ET AL. (2022, 2023) EVAPORATIVE EXTINCTION & DEPOSITION
# =============================================================================

"""
    renaudo_evaporation_parameter(T_c::Real, RH::Real) -> Real

Evaporation parameter k [s/m^2] from Renaudo et al. (2022, 2023):
    k = 1.25e9 / (T_c - T_bh)
where T_bh is thermodynamic wet-bulb temperature calculated via Stull (2011).
"""
function renaudo_evaporation_parameter(T_c::Real, RH::Real)
    T_bh = stull_twb(T_c, RH)
    delta_T = max(0.05, float(T_c) - T_bh)
    return 1.25e9 / delta_T
end

"""
    renaudo_extinction_diameter(H::Real, T_c::Real, RH::Real; mu_g=1.81e-5, rho_d=997.0, rho_g=1.2, g=9.80665) -> Real

Critical initial cut-off diameter d_min [m] below which an evaporating droplet completely vanishes
in mid-air before landing from release height H [m] (Renaudo et al. 2022 Eq. 10):
    d_min = [ 54 * mu_g * H / (k * (rho_d - rho_g) * g) ]^(1/4)
"""
function renaudo_extinction_diameter(H::Real, T_c::Real, RH::Real;
                                     mu_g::Real=1.81e-5,
                                     rho_d::Real=997.0,
                                     rho_g::Real=1.2,
                                     g::Real=9.80665)
    k = renaudo_evaporation_parameter(T_c, RH)
    delta_rho = max(1.0, float(rho_d) - float(rho_g))
    arg = (54.0 * float(mu_g) * float(H)) / (k * delta_rho * float(g))
    return (max(0.0, arg))^(0.25)
end

"""
    renaudo_deposition_time(d0::Real, H::Real, T_c::Real, RH::Real; kwargs...) -> Real

Total in-flight sedimentation time t_dep [s] until soil deposition (Renaudo et al. 2022 Eq. 8).
Returns Inf if d0 <= d_min (droplet completely evaporates in flight).
"""
function renaudo_deposition_time(d0::Real, H::Real, T_c::Real, RH::Real;
                                 mu_g::Real=1.81e-5,
                                 rho_d::Real=997.0,
                                 rho_g::Real=1.2,
                                 g::Real=9.80665)
    d_min = renaudo_extinction_diameter(H, T_c, RH; mu_g=mu_g, rho_d=rho_d, rho_g=rho_g, g=g)
    if float(d0) <= d_min
        return Inf
    end
    k = renaudo_evaporation_parameter(T_c, RH)
    ratio = (d_min / float(d0))^4
    return k * (float(d0)^2) * (1.0 - (1.0 - ratio)^(1.0 / 3.0))
end

"""
    renaudo_landing_diameter(d0::Real, H::Real, T_c::Real, RH::Real; kwargs...) -> Real

Remaining droplet diameter [m] upon impacting the ground from height H (Renaudo et al. 2022 Eq. 9).
Returns 0.0 if d0 <= d_min.
"""
function renaudo_landing_diameter(d0::Real, H::Real, T_c::Real, RH::Real;
                                 mu_g::Real=1.81e-5,
                                 rho_d::Real=997.0,
                                 rho_g::Real=1.2,
                                 g::Real=9.80665)
    d_min = renaudo_extinction_diameter(H, T_c, RH; mu_g=mu_g, rho_d=rho_d, rho_g=rho_g, g=g)
    if float(d0) <= d_min
        return 0.0
    end
    ratio = (d_min / float(d0))^4
    return float(d0) * ((1.0 - ratio)^(1.0 / 6.0))
end

# =============================================================================
# CHIU GROUP EVAPORATION NUMBER (Sirignano 2010, Chiu 1982)
# =============================================================================

"""
    chiu_group_evaporation_number(N_tot::Real, R_d::Real, R_c::Real; Re=0.0, Sc=0.6, Le=1.0) -> Real

Compute Chiu's Group Evaporation Number G (Sirignano 2010, Eq. 3.125):
    G = 3 * (1 + 0.276 * Re^(1/2) * Sc^(1/3)) * Le * N_tot^(2/3) * (R_d / R_c)
where:
- G < 0.1: Single isolated droplet regime
- 0.1 <= G <= 100: Internal group combustion/evaporation regime
- G > 100: External group / cloud boundary sheath regime
"""
function chiu_group_evaporation_number(N_tot::Real, R_d::Real, R_c::Real;
                                       Re::Real=0.0,
                                       Sc::Real=0.6,
                                       Le::Real=1.0)
    if R_c <= 0.0 || N_tot <= 0.0
        return 0.0
    end
    corr = 1.0 + 0.276 * sqrt(max(0.0, float(Re))) * (float(Sc)^(1.0 / 3.0))
    return 3.0 * corr * float(Le) * (float(N_tot)^(2.0 / 3.0)) * (float(R_d) / float(R_c))
end

"""
    chiu_evaporation_regime(G::Real) -> Symbol

Classify spray cloud evaporation regime based on Chiu's Group number:
- `:single_droplet` (G < 0.1)
- `:internal_group` (0.1 <= G <= 100)
- `:external_group` (G > 100)
"""
function chiu_evaporation_regime(G::Real)
    g_val = float(G)
    if g_val < 0.1
        return :single_droplet
    elseif g_val <= 100.0
        return :internal_group
    else
        return :external_group
    end
end

end # module Correlations
