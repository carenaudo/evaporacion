"""
    algebraic.jl

Calibrated empirical and physics-informed surrogate coefficients (Modelos 2b, 4, 5, 6, 6b, 6b_refit, 7).
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

function calc_C_model6b_refit(d_mm::T, U::Real, T_inf::Real; coeffs=nothing) where {T<:Real}
    Re = calc_reynolds_surrogate(d_mm, U, T_inf)
    t_ratio = max(T(T_inf) / T(25.0), T(1e-6))
    if coeffs === nothing
        # Calibrated directly against Julia Rigorous V2 ODE lifetimes (2,520 cases)
        if Re < T(10.0)
            A0, n0, A1, n1, ci = T(5.528445e-05), T(0.272604), T(9.536017e-05), T(0.159309), T(0.495)
        elseif Re < T(100.0)
            A0, n0, A1, n1, ci = T(5.495803e-05), T(0.270657), T(9.686232e-05), T(0.160461), T(0.496)
        else
            A0, n0, A1, n1, ci = T(5.356565e-05), T(0.266251), T(9.814937e-05), T(0.164165), T(0.493)
        end
        return A0 * (t_ratio^n0) + A1 * (t_ratio^n1) * ((d_mm * abs(T(U)))^ci)
    else
        # coeffs is a NamedTuple (; low=(A0,n0,A1,n1,ci), mid=..., high=...)
        reg = Re < T(10.0) ? coeffs.low : (Re < T(100.0) ? coeffs.mid : coeffs.high)
        A0, n0, A1, n1, ci = T(reg[1]), T(reg[2]), T(reg[3]), T(reg[4]), T(reg[5])
        return A0 * (t_ratio^n0) + A1 * (t_ratio^n1) * ((d_mm * abs(T(U)))^ci)
    end
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

function calc_C_model7(d_mm::T, U::Real, T_inf::Real, HR::Real, T_wb::Real, T0::Real; coeffs=nothing) where {T<:Real}
    Re = calc_reynolds_surrogate(d_mm, U, T_inf)
    t_ratio = max(T(T_inf) / T(25.0), T(1e-6))
    dT = T(T_inf) - T(T_wb)
    Theta = abs(dT) > T(1e-3) ? (T(T0) - T(T_wb)) / dT : zero(T)

    if coeffs === nothing
        scale = one(T)
        if Re < T(10.0)
            A0, n0, alpha = T(5.346202e-05), T(0.281287), T(-0.050473)
            A1, n1, beta  = T(9.376442e-05), T(0.167952), T(-0.028636)
            gamma, delta, ci = T(0.008057), T(0.000816), T(0.495)
        elseif Re < T(100.0)
            A0, n0, alpha = T(5.292208e-05), T(0.281515), T(-0.054463)
            A1, n1, beta  = T(9.533896e-05), T(0.166045), T(-0.027342)
            gamma, delta, ci = T(0.007001), T(0.000643), T(0.496)
        else
            A0, n0, alpha = T(5.145208e-05), T(0.283085), T(-0.056330)
            A1, n1, beta  = T(9.684243e-05), T(0.167299), T(-0.023637)
            gamma, delta, ci = T(0.006643), T(0.000690), T(0.493)
        end
    else
        scale = one(T)
        reg = Re < T(10.0) ? coeffs.low : (Re < T(100.0) ? coeffs.mid : coeffs.high)
        A0, n0, alpha, A1, n1, beta, gamma, delta, ci = T.(reg)
    end

    hr_clamped = clamp(T(HR), zero(T), T(0.999))
    f_HR = (one(T) - hr_clamped)^alpha
    g_HR = (one(T) - hr_clamped)^beta
    h_Theta = exp(gamma * Theta + delta * (Theta^2))

    base_C = A0 * (t_ratio^n0) * f_HR + A1 * (t_ratio^n1) * ((d_mm * abs(T(U)))^ci) * g_HR
    return scale * base_C * h_Theta
end
