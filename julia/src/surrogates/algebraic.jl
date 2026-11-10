"""
    algebraic.jl

Calibrated empirical and physics-informed surrogate coefficients (Modelos 2b, 4, 5, 6, 6b, 6b_refit, 7, 8, 8_dynamic).
All equations are generic over `T <: Real`.
"""

# =============================================================================
# REYNOLDS NUMBER FOR SURROGATES
# =============================================================================

function calc_reynolds_surrogate(d_mm::T, U::Real, T_inf::Real) where {T<:Real}
    d_m = d_mm / T(1000.0)
    T_inf_K = T(T_inf) + T(T_ZERO_C)
    rho = gas_density_air(T_inf_K)
    mu  = gas_viscosity_air(T_inf_K)
    return rho * abs(T(U)) * d_m / max(mu, T(1e-12))
end

# =============================================================================
# EMPIRICAL VOLUME MODELS 4, 5, 6, 6B, 6B_REFIT
# =============================================================================

function calc_C_model4(d_mm::T, U::Real, T_inf::Real) where {T<:Real}
    Re = calc_reynolds_surrogate(d_mm, U, T_inf)
    c = Re < T(10.0) ? T(0.834) : (Re < T(100.0) ? T(0.646) : T(0.563))
    A = T(4.8775e-5)
    T_inf_K = T(T_inf) + T(T_ZERO_C)
    D_T = gas_diffusivity_water_air(T_inf_K)
    D_25 = gas_diffusivity_water_air(T(298.15))
    return A * (D_T / D_25) * (one(T) + (d_mm * abs(T(U)))^c)
end

function calc_C_model5(d_mm::T, U::Real, T_inf::Real) where {T<:Real}
    Re = calc_reynolds_surrogate(d_mm, U, T_inf)
    if Re < T(10.0)
        A_i, c_i = T(4.267e-5), T(0.4849)
    elseif Re < T(100.0)
        A_i, c_i = T(4.885e-5), T(0.6413)
    else
        A_i, c_i = T(4.998e-5), T(0.5467)
    end
    t_ratio = max(T(T_inf) / T(25.0), T(1e-6))
    return A_i * (t_ratio^T(0.1598)) * (one(T) + (d_mm * abs(T(U)))^c_i)
end

function calc_C_model6(d_mm::T, U::Real, T_inf::Real) where {T<:Real}
    Re = calc_reynolds_surrogate(d_mm, U, T_inf)
    c_i = Re < T(10.0) ? T(0.495) : (Re < T(100.0) ? T(0.496) : T(0.493))
    A0, n0 = T(3.7810e-5), T(0.2095)
    A1, n1 = T(6.0292e-5), T(0.1382)
    t_ratio = max(T(T_inf) / T(25.0), T(1e-6))
    return A0 * (t_ratio^n0) + A1 * (t_ratio^n1) * ((d_mm * abs(T(U)))^c_i)
end

function calc_C_model6b(d_mm::T, U::Real, T_inf::Real) where {T<:Real}
    Re = calc_reynolds_surrogate(d_mm, U, T_inf)
    t_ratio = max(T(T_inf) / T(25.0), T(1e-6))
    if Re < T(10.0)
        A0, n0 = T(6.282876e-05), T(0.093366)
        A1, n1 = T(2.381775e-05), T(0.733521)
        ci = T(0.495)
    elseif Re < T(100.0)
        A0, n0 = T(6.632398e-05), T(0.164756)
        A1, n1 = T(2.270989e-05), T(-0.257162)
        ci = T(0.496)
    else
        A0, n0 = T(7.425671e-05), T(0.027836)
        A1, n1 = T(1.729154e-05), T(-0.057516)
        ci = T(0.493)
    end
    return A0 * (t_ratio^n0) + A1 * (t_ratio^n1) * ((d_mm * abs(T(U)))^ci)
end

function calc_C_model6b_refit(d_mm::T, U::Real, T_inf::Real) where {T<:Real}
    Re = calc_reynolds_surrogate(d_mm, U, T_inf)
    t_ratio = max(T(T_inf) / T(25.0), T(1e-6))
    scale = T(1.6366)
    if Re < T(10.0)
        A0, n0, A1, n1, ci = T(3.6681e-05), T(0.2499), T(5.0043e-05), T(0.1384), T(0.495)
    elseif Re < T(100.0)
        A0, n0, A1, n1, ci = T(3.6917e-05), T(0.2456), T(4.9725e-05), T(0.1379), T(0.496)
    else
        A0, n0, A1, n1, ci = T(3.6789e-05), T(0.2331), T(4.9754e-05), T(0.1459), T(0.493)
    end
    return scale * (A0 * (t_ratio^n0) + A1 * (t_ratio^n1) * ((d_mm * abs(T(U)))^ci))
end

function calc_C_model2b(d_mm::T, U::Real, T_inf::Real, T_wb::Real) where {T<:Real}
    T_avg = (T(T_inf) + T(T_wb)) / T(2.0)
    T_avg_K = T_avg + T(T_ZERO_C)
    T_wb_K  = T(T_wb) + T(T_ZERO_C)

    d_m = d_mm / T(1000.0)
    rho = gas_density_air(T_avg_K)
    mu  = gas_viscosity_air(T_avg_K)
    k_g = gas_thermal_conductivity_air(T_avg_K)
    cp  = gas_cp_air(T_avg_K)
    Pr  = mu * cp / k_g
    Re  = rho * abs(T(U)) * d_m / max(mu, T(1e-12))
    Nu  = T(2.0) + T(0.552) * sqrt(max(Re, zero(T))) * (Pr^T(1.0 / 3.0))

    rho_l = liquid_density_water(T_wb_K)
    h_v   = latent_heat_water(T_wb_K)

    C_m2_sK = sqrt(T(π)) * k_g * Nu / (rho_l * h_v)
    return C_m2_sK * T(1e6) # in mm^2/(s * K)
end

function calc_C_model7(d_mm::T, U::Real, T_inf::Real, HR::Real, T_wb::Real, T0::Real) where {T<:Real}
    Re = calc_reynolds_surrogate(d_mm, U, T_inf)
    t_ratio = max(T(T_inf) / T(25.0), T(1e-6))
    dT = T(T_inf) - T(T_wb)
    Theta = abs(dT) > T(1e-3) ? (T(T0) - T(T_wb)) / dT : zero(T)
    scale = T(1.6335)

    if Re < T(10.0)
        A0, n0, alpha = T(3.5677e-05), T(0.2615), T(-0.0427)
        A1, n1, beta  = T(4.9687e-05), T(0.1489), T(-0.0158)
        gamma, delta, ci = T(0.0071), T(0.0004), T(0.495)
    elseif Re < T(100.0)
        A0, n0, alpha = T(3.5902e-05), T(0.2576), T(-0.0420)
        A1, n1, beta  = T(4.9178e-05), T(0.1475), T(-0.0200)
        gamma, delta, ci = T(0.0061), T(0.0002), T(0.496)
    else
        A0, n0, alpha = T(3.5611e-05), T(0.2491), T(-0.0463)
        A1, n1, beta  = T(4.9417e-05), T(0.1529), T(-0.0147)
        gamma, delta, ci = T(0.0057), T(0.0002), T(0.493)
    end

    hr_clamped = clamp(T(HR), zero(T), T(0.999))
    f_HR = (one(T) - hr_clamped)^alpha
    g_HR = (one(T) - hr_clamped)^beta
    h_Theta = one(T) + gamma * Theta + delta * (Theta^2)

    base_C = A0 * (t_ratio^n0) * f_HR + A1 * (t_ratio^n1) * ((d_mm * abs(T(U)))^ci) * g_HR
    return scale * base_C * h_Theta
end

# =============================================================================
# COUPLED SURROGATES 8 & 8 DYNAMIC (dV/dt, dT/dt)
# =============================================================================

function calc_params_model8(d0_mm::T, U::Real, Tinf::Real, HR::Real, Twb::Real, T0::Real) where {T<:Real}
    Re0 = calc_reynolds_surrogate(d0_mm, U, Tinf)
    dT_pot = T(Tinf) - T(Twb)
    if abs(dT_pot) < T(1e-3)
        dT_pot = T(1e-3)
    end
    Theta0 = (T(T0) - T(Twb)) / dT_pot

    K1_base = T(6.367e-03) * T(0.9202)
    hr_c = clamp(T(HR), zero(T), T(0.999))
    K1 = K1_base * (one(T) + T(0.019) * sqrt(max(Re0, zero(T)))) * ((one(T) - hr_c)^T(0.734)) * (one(T) - T(0.131) * Theta0)

    K2_term = one(T) - T(19.975) * Theta0
    K2_term = clamp(K2_term, T(0.01), T(10.0))
    K2 = T(1.597e-03) * (one(T) - T(0.027) * sqrt(max(Re0, zero(T)))) * ((one(T) - hr_c)^T(-0.301)) * K2_term

    lam = T(20.0)
    CV = T(0.5)
    CT = T(-0.5)

    return (K1=K1, K2=K2, lam=lam, CV=CV, CT=CT)
end

function calc_params_model8_dynamic(
    V::T, Temp::Real, U::Real, Tinf::Real, HR::Real, Twb::Real;
    params=nothing
) where {T<:Real}
    d_mm = sphere_diameter_from_volume(V)
    Re = calc_reynolds_surrogate(d_mm, U, Tinf)
    dT_pot = T(Tinf) - T(Twb)
    if abs(dT_pot) < T(1e-3)
        dT_pot = T(1e-3)
    end
    Theta = (T(Temp) - T(Twb)) / dT_pot

    if params === nothing
        A1, b1, alpha1, beta1, CV = T(5.883990e-03), T(1.412334e-02), T(7.066836e-01), T(-9.941564e-02), T(4.365568e-01)
        A2, b2, alpha2, beta2, CT = T(2.114610e-02), T(3.731696e-02), T(-3.036976e-01), T(-2.000011e+01), T(-5.076667e-01)
        lam_0, gamma = T(20.0), T(1.760136e-04)
    else
        A1, b1, alpha1, beta1, CV, A2, b2, alpha2, beta2, CT, lam_0, gamma = T.(params)
    end

    hr_c = clamp(T(HR), zero(T), T(0.999))
    K1 = A1 * (one(T) + b1 * sqrt(max(Re, zero(T)))) * ((one(T) - hr_c)^alpha1) * (one(T) + beta1 * Theta)

    theta_term2 = clamp(one(T) + beta2 * Theta, T(0.01), T(10.0))
    K2 = A2 * (one(T) + b2 * sqrt(max(Re, zero(T)))) * ((one(T) - hr_c)^alpha2) * theta_term2

    lam = lam_0 * (one(T) + gamma * Theta)

    return (K1=K1, K2=K2, lam=lam, CV=CV, CT=CT)
end
