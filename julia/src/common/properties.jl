"""
    Evaporacion.Common.Properties

High-Fidelity (Level 3) thermophysical properties for water and air.
All functions are generic on `T <: Real` to support automatic differentiation
(e.g., ForwardDiff.Dual) and high-performance compilation.
All internal calculations use SI units (K, Pa, kg, m, s, J, W).
"""
module Properties

using ..Constants: R_UNIVERSAL, M_AIR, M_WATER, P_ATM, T_ZERO_C, T_CRIT_WATER

export psat_water, liquid_density_water, liquid_cp_water, latent_heat_water,
       liquid_viscosity_water, liquid_thermal_conductivity_water,
       liquid_thermal_expansion_water,
       gas_viscosity_air, gas_thermal_conductivity_air, gas_diffusivity_water_air,
       gas_density_air, gas_cp_air, vapor_cp_water,
       vapor_viscosity_water, vapor_thermal_conductivity_water,
       wilke_mixture_viscosity, wilke_mixture_conductivity,
       reference_state_one_third, mole_fraction_from_mass_fraction, mass_fraction_from_mole_fraction,
       gas_mixture_properties,
       surface_mass_fraction_from_psat, stull_twb

# =============================================================================
# SATURATION VAPOR PRESSURE OF WATER [Pa]
# =============================================================================

"""
    psat_water(T_K::Real) -> Real

Saturation vapor pressure of liquid water in Pa using the Buck (1981) equation.
Accurate over 273.15 K to 373.15 K.
"""
function psat_water(T_K::T) where {T<:Real}
    T_C = T_K - T(T_ZERO_C)
    # Buck (1981) formula
    num = (T(18.678) - T_C / T(234.5)) * T_C
    den = T(257.14) + T_C
    return T(611.21) * exp(num / den)
end

# =============================================================================
# LIQUID WATER PROPERTIES
# =============================================================================

"""
    liquid_density_water(T_K::Real) -> Real

Liquid water density in kg/m^3 using Kell (1975) 5th-order polynomial.
"""
function liquid_density_water(T_K::T) where {T<:Real}
    t = T_K - T(T_ZERO_C) # temperature in Celsius
    return T(999.842594) + T(6.793952e-2) * t - T(9.095290e-3) * t^2 +
           T(1.001685e-4) * t^3 - T(1.120083e-6) * t^4 + T(6.536332e-9) * t^5
end

"""
    liquid_thermal_expansion_water(T_K::Real) -> Real

Volumetric thermal expansion coefficient of liquid water in 1/K:
    β = -(1/ρ) * (∂ρ/∂T)
computed from the derivative of Kell (1975).
"""
function liquid_thermal_expansion_water(T_K::T) where {T<:Real}
    t = T_K - T(T_ZERO_C)
    rho = liquid_density_water(T_K)
    drho_dT = T(6.793952e-2) - T(2.0) * T(9.095290e-3) * t +
              T(3.0) * T(1.001685e-4) * t^2 - T(4.0) * T(1.120083e-6) * t^3 +
              T(5.0) * T(6.536332e-9) * t^4
    return -drho_dT / rho
end

"""
    liquid_cp_water(T_K::Real) -> Real

Liquid water isobaric specific heat capacity in J/(kg * K).
"""
function liquid_cp_water(T_K::T) where {T<:Real}
    t = T_K - T(T_ZERO_C)
    return T(4217.4) - T(2.0) * t + T(0.01) * t^2
end

"""
    latent_heat_water(T_K::Real) -> Real

Enthalpy of vaporization (latent heat) of water in J/kg using the Watson correlation.
Reference state: T_0 = 373.15 K, ΔH_0 = 2.257e6 J/kg, T_crit = 647.096 K.
"""
function latent_heat_water(T_K::T) where {T<:Real}
    T_crit = T(T_CRIT_WATER)
    T_boil = T(373.15)
    # Ensure temperature does not exceed critical
    T_safe = min(T_K, T_crit - T(1e-4))
    tau = max((T_crit - T_safe) / (T_crit - T_boil), T(0.0))
    return T(2.257e6) * (tau^T(0.38))
end

"""
    liquid_viscosity_water(T_K::Real) -> Real

Dynamic viscosity of liquid water in Pa*s using Andrade / DIPPR equation.
"""
function liquid_viscosity_water(T_K::T) where {T<:Real}
    T_C = T_K - T(T_ZERO_C)
    return T(2.414e-5) * (T(10.0)^(T(247.8) / (T_C + T(133.15))))
end

"""
    liquid_thermal_conductivity_water(T_K::Real) -> Real

Thermal conductivity of liquid water in W/(m * K).
"""
function liquid_thermal_conductivity_water(T_K::T) where {T<:Real}
    T_C = T_K - T(T_ZERO_C)
    val = T(0.6065) - T(0.00122) * max(T_C - T(20.0), T(-20.0))
    return max(val, T(0.1))
end

# =============================================================================
# GAS (AIR & WATER VAPOR) PROPERTIES
# =============================================================================

"""
    gas_viscosity_air(T_K::Real) -> Real

Dynamic viscosity of air in Pa*s using Sutherland's law.
S0 = 110.4 K, T0 = 273.15 K, mu0 = 1.716e-5 Pa*s.
"""
function gas_viscosity_air(T_K::T) where {T<:Real}
    T_ref = T(273.15)
    S0    = T(110.4)
    mu0   = T(1.716e-5)
    return mu0 * ((T_K / T_ref)^T(1.5)) * (T_ref + S0) / (T_K + S0)
end

"""
    gas_thermal_conductivity_air(T_K::Real) -> Real

Thermal conductivity of air in W/(m * K).
"""
function gas_thermal_conductivity_air(T_K::T) where {T<:Real}
    return T(2.411e-2) + T(7.704e-5) * T_K - T(2.363e-8) * (T_K^2)
end

"""
    gas_diffusivity_water_air(T_K::Real, P_Pa::Real=P_ATM) -> Real

Binary diffusion coefficient of water vapor in air in m^2/s using Fuller, Schettler & Giddings (1966).
"""
function gas_diffusivity_water_air(T_K::T, P_Pa::Real=P_ATM) where {T<:Real}
    M_a = T(28.97)
    M_w = T(18.015)
    P_val = T(P_Pa)
    num = T(1.013e-2) * (T_K^T(1.75)) * sqrt(T(1.0) / M_a + T(1.0) / M_w)
    den = P_val * ((T(19.7)^T(1.0 / 3.0) + T(13.1)^T(1.0 / 3.0))^2)
    return num / den
end

"""
    gas_density_air(T_K::Real, P_Pa::Real=P_ATM) -> Real

Dry air density in kg/m^3 from ideal gas law.
"""
function gas_density_air(T_K::T, P_Pa::Real=P_ATM) where {T<:Real}
    return T(P_Pa) * T(M_AIR) / (T(R_UNIVERSAL) * T_K)
end

"""
    gas_cp_air(T_K::Real) -> Real

Isobaric specific heat capacity of dry air in J/(kg * K).
"""
function gas_cp_air(T_K::T) where {T<:Real}
    return T(1006.0)
end

"""
    vapor_cp_water(T_K::Real) -> Real

Isobaric specific heat capacity of water vapor in J/(kg * K).
"""
function vapor_cp_water(T_K::T) where {T<:Real}
    return T(1880.0)
end

"""
    vapor_viscosity_water(T_K::Real) -> Real

Dynamic viscosity of low-pressure water vapor in Pa*s using Sutherland's law for steam:
    mu = mu0 * (T / T0)^(3/2) * (T0 + S) / (T + S)
with T0 = 273.15 K, S = 647.0 K, mu0 = 8.85e-6 Pa*s.
"""
function vapor_viscosity_water(T_K::T) where {T<:Real}
    T_ref = T(273.15)
    S0    = T(647.0)
    mu0   = T(8.85e-6)
    return mu0 * ((T_K / T_ref)^T(1.5)) * (T_ref + S0) / (T_K + S0)
end

"""
    vapor_thermal_conductivity_water(T_K::Real) -> Real

Thermal conductivity of low-pressure water vapor in W/(m * K).
"""
function vapor_thermal_conductivity_water(T_K::T) where {T<:Real}
    T_C = T_K - T(T_ZERO_C)
    return T(0.0176) + T(5.56e-5) * T_C
end

"""
    reference_state_one_third(T_s::Real, T_inf::Real, Y_v_s::Real, Y_v_inf::Real) -> (T_ref, Y_v_ref)

Computes the 1/3 reference state (Hubbard, Denny & Mills, 1975; Abramzon & Sirignano, 1989):
    T_ref = T_s + (1/3) * (T_inf - T_s)
    Y_v_ref = Y_v_s + (1/3) * (Y_v_inf - Y_v_s)
"""
function reference_state_one_third(T_s::T, T_inf::Real, Y_v_s::Real, Y_v_inf::Real) where {T<:Real}
    one_third = T(1.0 / 3.0)
    T_ref   = T_s + one_third * (T(T_inf) - T_s)
    Y_v_ref = T(Y_v_s) + one_third * (T(Y_v_inf) - T(Y_v_s))
    return (T_ref = T_ref, Y_v_ref = Y_v_ref)
end

"""
    mole_fraction_from_mass_fraction(Y_v::Real; M_v::Real=M_WATER, M_g::Real=M_AIR) -> Real

Convert vapor mass fraction to mole fraction in a binary vapor-air mixture.
"""
function mole_fraction_from_mass_fraction(Y_v::T; M_v::Real=M_WATER, M_g::Real=M_AIR) where {T<:Real}
    y = clamp(Y_v, zero(T), T(0.99999999))
    num = y / T(M_v)
    den = num + (one(T) - y) / T(M_g)
    return num / den
end

"""
    mass_fraction_from_mole_fraction(X_v::Real; M_v::Real=M_WATER, M_g::Real=M_AIR) -> Real

Convert vapor mole fraction to mass fraction in a binary vapor-air mixture.
"""
function mass_fraction_from_mole_fraction(X_v::T; M_v::Real=M_WATER, M_g::Real=M_AIR) where {T<:Real}
    x = clamp(X_v, zero(T), T(0.99999999))
    num = x * T(M_v)
    den = num + (one(T) - x) * T(M_g)
    return num / den
end

"""
    wilke_mixture_viscosity(mu_v::Real, mu_a::Real, X_v::Real, X_a::Real; M_v::Real=M_WATER, M_a::Real=M_AIR) -> Real

Dynamic viscosity of a binary gas mixture using Wilke's semi-empirical mixing rule (1950):
    mu_mix = X_v * mu_v / (X_v + X_a * Phi_va) + X_a * mu_a / (X_a + X_v * Phi_av)
"""
function wilke_mixture_viscosity(
    mu_v::T, mu_a::Real, X_v::Real, X_a::Real;
    M_v::Real=M_WATER, M_a::Real=M_AIR
) where {T<:Real}
    mv = T(M_v)
    ma = T(M_a)
    mu_1 = mu_v
    mu_2 = T(mu_a)
    x_1  = T(X_v)
    x_2  = T(X_a)

    phi_12 = ((one(T) + sqrt(mu_1 / mu_2) * ((ma / mv)^T(0.25)))^2) / (sqrt(T(8.0) * (one(T) + mv / ma)))
    phi_21 = ((one(T) + sqrt(mu_2 / mu_1) * ((mv / ma)^T(0.25)))^2) / (sqrt(T(8.0) * (one(T) + ma / mv)))

    denom1 = x_1 + x_2 * phi_12
    denom2 = x_2 + x_1 * phi_21
    term1 = denom1 > zero(T) ? (x_1 * mu_1 / denom1) : zero(T)
    term2 = denom2 > zero(T) ? (x_2 * mu_2 / denom2) : zero(T)
    return term1 + term2
end

"""
    wilke_mixture_conductivity(k_v::Real, k_a::Real, X_v::Real, X_a::Real; mu_v::Real=0.0, mu_a::Real=0.0, M_v::Real=M_WATER, M_a::Real=M_AIR) -> Real

Thermal conductivity of a binary gas mixture using the Wilke-Mason-Saxena formula
with interaction parameters Phi_12 and Phi_21 computed from viscosities.
"""
function wilke_mixture_conductivity(
    k_v::T, k_a::Real, X_v::Real, X_a::Real;
    mu_v::Real=0.0, mu_a::Real=0.0, M_v::Real=M_WATER, M_a::Real=M_AIR
) where {T<:Real}
    mv = T(M_v)
    ma = T(M_a)
    k_1 = k_v
    k_2 = T(k_a)
    x_1 = T(X_v)
    x_2 = T(X_a)

    # Use viscosity ratio for Phi_ij if provided, else conductivity ratio
    ratio = (mu_v > 0.0 && mu_a > 0.0) ? (T(mu_v) / T(mu_a)) : (k_1 / k_2)
    phi_12 = ((one(T) + sqrt(ratio) * ((ma / mv)^T(0.25)))^2) / (sqrt(T(8.0) * (one(T) + mv / ma)))
    phi_21 = ((one(T) + sqrt(one(T) / ratio) * ((mv / ma)^T(0.25)))^2) / (sqrt(T(8.0) * (one(T) + ma / mv)))

    denom1 = x_1 + x_2 * phi_12
    denom2 = x_2 + x_1 * phi_21
    term1 = denom1 > zero(T) ? (x_1 * k_1 / denom1) : zero(T)
    term2 = denom2 > zero(T) ? (x_2 * k_2 / denom2) : zero(T)
    return term1 + term2
end

"""
    gas_mixture_properties(T_ref_K::Real, Y_v_ref::Real, P_Pa::Real=P_ATM)

Calculates the complete thermophysical properties of a humid air / water vapor mixture
evaluated at reference temperature `T_ref_K` and reference vapor mass fraction `Y_v_ref`
using the Wilke binary formulation.
Returns `(rho_g, mu_g, k_g, D_ab, cp_g, M_mix)`.
"""
function gas_mixture_properties(T_ref_K::T, Y_v_ref::Real, P_Pa::Real=P_ATM) where {T<:Real}
    P_val = T(P_Pa)
    y_v = clamp(T(Y_v_ref), zero(T), T(0.99999999))
    y_a = one(T) - y_v

    # Mole fractions
    x_v = mole_fraction_from_mass_fraction(y_v)
    x_a = one(T) - x_v

    # Apparent molar mass
    M_mix = x_v * T(M_WATER) + x_a * T(M_AIR)

    # Ideal gas density of mixture
    rho_g = P_val * M_mix / (T(R_UNIVERSAL) * T_ref_K)

    # Pure species properties
    mu_a = gas_viscosity_air(T_ref_K)
    mu_v = vapor_viscosity_water(T_ref_K)
    k_a  = gas_thermal_conductivity_air(T_ref_K)
    k_v  = vapor_thermal_conductivity_water(T_ref_K)
    cp_a = gas_cp_air(T_ref_K)
    cp_v = vapor_cp_water(T_ref_K)

    # Wilke mixing for viscosity and conductivity
    mu_g = wilke_mixture_viscosity(mu_v, mu_a, x_v, x_a)
    k_g  = wilke_mixture_conductivity(k_v, k_a, x_v, x_a; mu_v=mu_v, mu_a=mu_a)

    # Mixture specific heat (mass-weighted)
    cp_g = y_v * cp_v + y_a * cp_a

    # Binary diffusion coefficient
    D_ab = gas_diffusivity_water_air(T_ref_K, P_val)

    return (
        rho_g = rho_g,
        mu_g  = mu_g,
        k_g   = k_g,
        D_ab  = D_ab,
        cp_g  = cp_g,
        M_mix = M_mix
    )
end

# =============================================================================
# CONVERSIONS & PSYCHROMETRICS
# =============================================================================

"""
    surface_mass_fraction_from_psat(P_sat_Pa::Real, P_total_Pa::Real=P_ATM) -> Real

Calculate mass fraction of water vapor at the liquid surface from saturation vapor pressure.
"""
function surface_mass_fraction_from_psat(P_sat_Pa::T, P_total_Pa::Real=P_ATM) where {T<:Real}
    P_tot = T(P_total_Pa)
    P_vap = clamp(P_sat_Pa, T(0.0), T(0.999999) * P_tot)
    denom = P_tot - P_vap
    if denom <= T(0.0)
        return T(0.999999999)
    end
    ratio = (T(M_WATER) / T(M_AIR)) * P_vap / denom
    return ratio / (T(1.0) + ratio)
end

"""
    stull_twb(T_C::Real, RH_pct::Real) -> Real

Stull (2011) empirical wet-bulb temperature formula in Celsius.
Valid for relative humidity in [1, 99]% and ambient temperature in [-20, 50] °C.
"""
function stull_twb(T_C::T, RH_pct::Real) where {T<:Real}
    rh = T(RH_pct)
    t = T_C
    return (
        t * atan(T(0.151977) * sqrt(rh + T(8.313659)))
        + atan(t + rh)
        - atan(rh - T(1.676331))
        + T(0.00391838) * (rh^T(1.5)) * atan(T(0.023101) * rh)
        - T(4.686035)
    )
end

end # module Properties
