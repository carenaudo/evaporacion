"""
Surface Area Analysis Across All Cases & Models
Auxiliary Material Script

Calculates the time-averaged mean surface area:
  S_mean = (1 / t_evap) * ∫_0^{t_evap} S(t) dt
and the temporal median surface area:
  S_median = S(t = t_evap / 2)
for:
  1. In-flight droplets (all 9 cases × 4 diameters × 10 models)
  2. Sessile droplets (all 9 cases × 3 contact angles × 8 models)
  3. Cross-domain comparison for d0 = 200 µm

Exports CSV tables and appends to the auxiliary material datasets.
"""

ENV["GKSwstype"] = "100"
using DataFrames
using CSV
using Printf
using QuadGK

# Load Evaporacion.jl package if not already loaded
if !isdefined(Main, :Evaporacion)
    include(joinpath(@__DIR__, "..", "..", "..", "julia", "src", "Evaporacion.jl"))
    using .Evaporacion
end

const SCRIPT_DIR = @__DIR__
const DATA_DIR = normpath(joinpath(SCRIPT_DIR, "..", "data"))
mkpath(DATA_DIR)

println("================================================================================")
println("COMPUTING MEAN & MEDIAN SURFACE AREAS ACROSS ALL CASES & MODELS")
println("================================================================================")

# 1. Load Environmental Cases
const CLICAP_CASES = [
    (id = 1, T_K = 283.15, T_C = 10.0, HR = 20.0),
    (id = 2, T_K = 283.15, T_C = 10.0, HR = 50.0),
    (id = 3, T_K = 283.15, T_C = 10.0, HR = 80.0),
    (id = 4, T_K = 298.15, T_C = 25.0, HR = 20.0),
    (id = 5, T_K = 298.15, T_C = 25.0, HR = 50.0),
    (id = 6, T_K = 298.15, T_C = 25.0, HR = 80.0),
    (id = 7, T_K = 313.15, T_C = 40.0, HR = 20.0),
    (id = 8, T_K = 313.15, T_C = 40.0, HR = 50.0),
    (id = 9, T_K = 313.15, T_C = 40.0, HR = 80.0)
]

const DIAMETERS_UM = [50.0, 200.0, 600.0, 1200.0]
const P_ATM = 101325.0

# -----------------------------------------------------------------------------
# HELPER: Compute Mean and Median from Discrete Time Series (t, S)
# -----------------------------------------------------------------------------
function calc_mean_and_median_area(t_vec::AbstractVector{<:Real}, S_vec::AbstractVector{<:Real})
    N = length(t_vec)
    if N < 2
        return (S_vec[1], S_vec[1])
    end

    t_end = t_vec[end]
    t_start = t_vec[1]
    dt_total = t_end - t_start

    if dt_total <= 0.0
        return (S_vec[1], S_vec[1])
    end

    # Trapezoidal integration for time-averaged mean
    integral = 0.0
    for i in 1:(N - 1)
        dt = t_vec[i+1] - t_vec[i]
        integral += 0.5 * (S_vec[i] + S_vec[i+1]) * dt
    end
    S_mean = integral / dt_total

    # Median over time = value at t = t_start + dt_total / 2
    t_half = t_start + 0.5 * dt_total
    # Linear interpolation at t_half
    idx = searchsortedlast(t_vec, t_half)
    idx = clamp(idx, 1, N - 1)
    frac = (t_half - t_vec[idx]) / max(t_vec[idx+1] - t_vec[idx], 1e-12)
    S_median = S_vec[idx] + frac * (S_vec[idx+1] - S_vec[idx])

    return (S_mean, S_median)
end

# -----------------------------------------------------------------------------
# PART 1: IN-FLIGHT DROPLETS
# -----------------------------------------------------------------------------
println("\n[1/2] Computing surface areas for in-flight droplets...")

records_inf_sa = []

for case in CLICAP_CASES
    T_C = case.T_C
    T_K = case.T_K
    HR_pct = case.HR
    HR_frac = HR_pct / 100.0
    psat_inf = psat_water(T_K)
    pv_inf = HR_frac * psat_inf
    M_w = Evaporacion.Common.Constants.M_WATER
    M_a = Evaporacion.Common.Constants.M_AIR
    Y_inf = (pv_inf * M_w) / (pv_inf * M_w + (P_ATM - pv_inf) * M_a)
    T_wb_C = stull_twb(T_C, HR_pct)

    for d_um in DIAMETERS_UM
        d_m = d_um * 1e-6
        S0_mm2 = π * (d_um * 1e-3)^2

        for (m_key, m_id, m_prof, m_desc) in [
            ("a", :d2_classic, :uniforme, "a (D² Clásica, Sh=2)"),
            ("b", :d2_stefan, :uniforme, "b (D² con Stefan)"),
            ("c", :d2_convective, :uniforme, "c (D² con Convección)"),
            ("d", :abramzon_sirignano, :parabolico, "d (A-S Parabólico)"),
            ("d_uni", :abramzon_sirignano, :uniforme, "d_uni (A-S Uniforme)"),
            ("d_zhou", :abramzon_sirignano, :zhou_cubic, "d_zhou (A-S Cúbico Zhou)"),
            ("revised", :d2_revised_dallabarba, :uniforme, "revised (D² Dalla Barba c/ Twb ventilado)"),
            ("vlasov", :vlasov_kinetic, :uniforme, "vlasov (Vlasov 2021 Cinético-Difusivo)"),
            ("surr_6b", :surrogate_6b, :uniforme, "surr_6b (Surrogante 6b Refit)"),
            ("surr_8dyn", :surrogate_8dyn, :uniforme, "surr_8dyn (Surrogante 8 Dinámico)")
        ]
            t_evap = 0.0
            S_mean = 0.0
            S_med  = 0.0

            if m_id in (:d2_classic, :d2_stefan, :d2_convective, :d2_revised_dallabarba, :abramzon_sirignano)
                isothermal = m_id in (:d2_classic, :d2_stefan, :d2_convective, :d2_revised_dallabarba)
                is_conv = m_id in (:d2_convective, :d2_revised_dallabarba, :abramzon_sirignano)
                T_init_C = isothermal ? T_wb_C : T_C

                rho_l0 = liquid_density_water(T_init_C + 273.15)
                m0 = rho_l0 * (π / 6.0) * (d_m^3)

                params = DropletParams(
                    T_inf = T_C,
                    Y_inf = Y_inf,
                    P_atm = P_ATM,
                    model_variant = m_id,
                    profile = m_prof,
                    terminal_velocity = is_conv,
                    isothermal = isothermal,
                    T_d_fixed = T_wb_C
                )

                t_est = (d_um / 50.0)^2 * 10.0
                tspan = (0.0, max(50.0, t_est * 4.0))

                sol = solve_droplet(params, [m0, T_init_C], tspan)
                t_evap = sol.t[end]

                # Reconstruct S(t) = π * d(t)^2 from sol
                S_t = Float64[]
                for state in sol.u
                    m_curr = max(state[1], 1e-18)
                    T_curr = state[2]
                    rho_curr = liquid_density_water(T_curr + 273.15)
                    d_curr = cbrt(m_curr / (rho_curr * π / 6.0))
                    push!(S_t, π * (d_curr * 1e3)^2) # in mm^2
                end

                S_mean, S_med = calc_mean_and_median_area(sol.t, S_t)

            elseif m_id == :vlasov_kinetic
                # Vlasov analytical trajectory R(t) = -l_k + sqrt((R0+l_k)^2 - 2*D*M/rho * delta_c * t)
                rho_l0 = liquid_density_water(T_wb_C + 273.15)
                D_val = gas_diffusivity_water_air(T_K, P_ATM)
                delta_c = (psat_water(T_wb_C + 273.15) / (Evaporacion.Common.Constants.R_GAS * (T_wb_C + 273.15)) * M_w) * (1.0 - HR_frac)
                R0 = d_m / 2.0
                alpha_evap = 0.04
                v_th = sqrt(8.0 * Evaporacion.Common.Constants.R_GAS * T_K / (π * M_w))

                tau_diff = (rho_l0 * R0^2) / (2.0 * D_val * delta_c)
                l_k = vlasov_kinetic_length(D_val, alpha_evap, v_th)
                t_evap = tau_diff * (1.0 + 2.0 * l_k / R0)

                t_grid = range(0.0, t_evap, length=100)
                S_t = [4.0 * π * (vlasov_radius_trajectory(t, R0, D_val, delta_c, rho_l0, M_w, alpha_evap, v_th) * 1e3)^2 for t in t_grid]
                S_mean, S_med = calc_mean_and_median_area(collect(t_grid), S_t)

            elseif m_id in (:surrogate_6b, :surrogate_8dyn)
                stype = m_id == :surrogate_6b ? :modelo_6b_refit : :modelo_8_dynamic
                v_term_est = calc_terminal_velocity(d_m, T_K, P_ATM)
                cfg = SurrogateConfig(
                    model_type = stype,
                    d0_mm = d_um * 1e-3,
                    T_inf = T_C,
                    T_wb = T_wb_C,
                    U = max(0.2, v_term_est),
                    HR = HR_frac,
                    T0 = T_C,
                    t_end = 2500.0
                )
                sol = solve_surrogate(cfg)
                t_evap = sol.t[end]

                # V(t) in mm^3 -> S(t) = π * (6V/π)^(2/3)
                S_t = [π * ((6.0 * max(u[1], 1e-12) / π)^(2.0/3.0)) for u in sol.u]
                S_mean, S_med = calc_mean_and_median_area(sol.t, S_t)
            end

            push!(records_inf_sa, (
                Dominio = "En Vuelo",
                Caso = case.id,
                T_C = T_C,
                HR_pct = HR_pct,
                d0_um = d_um,
                Modelo_Key = m_key,
                Modelo_Desc = m_desc,
                t_evap_s = round(t_evap, digits=3),
                S0_mm2 = round(S0_mm2, digits=6),
                S_mean_mm2 = round(S_mean, digits=6),
                S_median_mm2 = round(S_med, digits=6),
                Ratio_Smean_S0 = round(S_mean / S0_mm2, digits=4),
                Ratio_Smed_S0 = round(S_med / S0_mm2, digits=4)
            ))
        end
    end
end

df_inf_sa = DataFrame(records_inf_sa)
CSV.write(joinpath(DATA_DIR, "surface_areas_inflight_all_cases.csv"), df_inf_sa)
println("Saved in-flight surface areas: $(nrow(df_inf_sa)) rows")

# -----------------------------------------------------------------------------
# PART 2: SESSILE DROPLETS
# -----------------------------------------------------------------------------
println("\n[2/2] Computing surface areas for sessile droplets...")

const CONTACT_ANGLES_DEG = [20.0, 40.0, 80.0]
const D0_SESSILE_UM = 200.0
const D0_M_SES = D0_SESSILE_UM * 1e-6
const V0_M3_SES = (π / 6.0) * (D0_M_SES^3)

records_ses_sa = []

for case in CLICAP_CASES
    T_C = case.T_C
    T_K = case.T_K
    HR_pct = case.HR
    HR_frac = HR_pct / 100.0

    rho_l = liquid_density_water(T_K)
    D_val = gas_diffusivity_water_air(T_K, P_ATM)
    psat = psat_water(T_K)
    c_sat = (psat * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * T_K)
    delta_c = (1.0 - HR_frac) * c_sat

    for th_deg in CONTACT_ANGLES_DEG
        th0_rad = deg2rad(th_deg)
        R0 = r_from_volume_sessile(V0_M3_SES, th0_rad)

        # Initial liquid-gas spherical cap surface area S_cap,0 = 2π R0^2 / (1 + cos(θ₀))
        S_cap0_mm2 = (2.0 * π * (R0 * 1e3)^2 / (1.0 + cos(th0_rad)))
        S_base0_mm2 = π * (R0 * 1e3)^2

        # 1. Mode AC (Constant Angle)
        # S_cap(t) = 2π R(t)^2 / (1 + cos(θ₀)) = S_cap,0 * (1 - t/t_evap)
        cfg_ca = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CA, t_end=2000.0)
        sol_ca = solve_sessile(cfg_ca)
        t_ca = sol_ca.t_evap_analytical
        # For pure linear R^2(t), S_mean = 0.5 * S_cap0, S_med = 0.5 * S_cap0
        S_mean_ca = 0.5 * S_cap0_mm2
        S_med_ca = 0.5 * S_cap0_mm2

        # 2. Mode RC (Constant Radius, Pinned)
        # R(t) = R0, θ(t) goes from θ₀ to 0. S_cap(t) = 2π R0^2 / (1 + cos(θ(t)))
        cfg_cr = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CR, t_end=2000.0)
        sol_cr = solve_sessile(cfg_cr)
        t_cr = sol_cr.t_evap_analytical

        # Reconstruct S_cap(t) from ODE trajectory or analytical relation
        S_cap_cr_t = [2.0 * π * (R0 * 1e3)^2 / (1.0 + cos(th)) for th in sol_cr.theta]
        S_mean_cr, S_med_cr = calc_mean_and_median_area(sol_cr.t, S_cap_cr_t)

        # 3. Mode Stick-Slide (SS)
        cfg_ss = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, theta_rec=th0_rad/2.0, mode=:SS, t_end=2000.0)
        sol_ss = solve_sessile(cfg_ss)
        t_ss = sol_ss.t_evap_analytical
        S_cap_ss_t = [2.0 * π * (sol_ss.R[i] * 1e3)^2 / (1.0 + cos(sol_ss.theta[i])) for i in 1:length(sol_ss.t)]
        S_mean_ss, S_med_ss = calc_mean_and_median_area(sol_ss.t, S_cap_ss_t)

        # 4. Mode Stick-Jump (SJ)
        cfg_sj = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, theta_min=th0_rad/2.0, theta_max=th0_rad, mode=:SJ, t_end=2000.0)
        sol_sj = solve_sessile(cfg_sj)
        t_sj = sol_sj.t_evap_analytical
        S_cap_sj_t = [2.0 * π * (sol_sj.R[i] * 1e3)^2 / (1.0 + cos(sol_sj.theta[i])) for i in 1:length(sol_sj.t)]
        S_mean_sj, S_med_sj = calc_mean_and_median_area(sol_sj.t, S_cap_sj_t)

        # 5. Substrate Cooled on Glass (CR mode)
        L_val = latent_heat_water(T_K)
        csat_plus = (psat_water(T_K + 0.5) * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * (T_K + 0.5))
        csat_minus = (psat_water(T_K - 0.5) * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * (T_K - 0.5))
        dcsat_dT_val = csat_plus - csat_minus
        tparams = SessileThermalParams(L=L_val, h_s=1e-3, k_s_eff=1.4, dcsat_dT=dcsat_dT_val)
        delta_c_eff = cooling_effective_delta_c(R0, th0_rad, D_val, delta_c, tparams)
        cfg_cool = SessileDropletConfig(D=D_val, delta_c=delta_c_eff, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CR, t_end=2000.0)
        sol_cool = solve_sessile(cfg_cool)
        t_cool = sol_cool.t_evap_analytical
        S_cap_cool_t = [2.0 * π * (R0 * 1e3)^2 / (1.0 + cos(th)) for th in sol_cool.theta]
        S_mean_cool, S_med_cool = calc_mean_and_median_area(sol_cool.t, S_cap_cool_t)

        # 6. Cooperative Vapor Shielding (Pair at B = 3 R0, CR mode)
        t_shielded = lifetime_pair_cr(3.0 * R0, R0, th0_rad, D_val, delta_c, rho_l)
        # Geometry matches CR mode scaled by shielding factor
        S_mean_shield = S_mean_cr
        S_med_shield  = S_med_cr

        for (m_lbl, t_val, s_mean, s_med) in [
            ("AC (CA)", t_ca, S_mean_ca, S_med_ca),
            ("RC (CR)", t_cr, S_mean_cr, S_med_cr),
            ("Stick-Slide (S-S)", t_ss, S_mean_ss, S_med_ss),
            ("Stick-Jump (S-J)", t_sj, S_mean_sj, S_med_sj),
            ("RC Substrate Cooled", t_cool, S_mean_cool, S_med_cool),
            ("RC Shielded (B=3R₀)", t_shielded, S_mean_shield, S_med_shield)
        ]
            push!(records_ses_sa, (
                Dominio = "Sésil",
                Caso = case.id,
                T_C = T_C,
                HR_pct = HR_pct,
                d0_equiv_um = D0_SESSILE_UM,
                theta0_deg = th_deg,
                R0_um = round(R0 * 1e6, digits=1),
                Modelo = m_lbl,
                t_evap_s = round(t_val, digits=3),
                S_cap0_mm2 = round(S_cap0_mm2, digits=6),
                S_base0_mm2 = round(S_base0_mm2, digits=6),
                S_mean_mm2 = round(s_mean, digits=6),
                S_median_mm2 = round(s_med, digits=6),
                Ratio_Smean_S0 = round(s_mean / S_cap0_mm2, digits=4),
                Ratio_Smed_S0 = round(s_med / S_cap0_mm2, digits=4)
            ))
        end
    end
end

df_ses_sa = DataFrame(records_ses_sa)
CSV.write(joinpath(DATA_DIR, "surface_areas_sessile_all_cases.csv"), df_ses_sa)
println("Saved sessile surface areas: $(nrow(df_ses_sa)) rows")

# -----------------------------------------------------------------------------
# PART 3: UNIFIED COMPARATIVE SUMMARY TABLE FOR PUBLICATION (d0 = 200 µm)
# -----------------------------------------------------------------------------
println("\n[3/3] Generating Master Publication Table of Mean & Median Surface Areas...")

records_summary = []

for c_id in 1:9
    c_info = CLICAP_CASES[c_id]
    t_c = c_info.T_C
    hr = c_info.HR

    # Representative In-Flight Models for d0 = 200 µm
    sub_inf = filter(r -> r.Caso == c_id && r.d0_um == 200.0, df_inf_sa)
    m_a  = filter(r -> r.Modelo_Key == "a", sub_inf)[1, :]
    m_c  = filter(r -> r.Modelo_Key == "c", sub_inf)[1, :]
    m_d  = filter(r -> r.Modelo_Key == "d", sub_inf)[1, :]
    m_zh = filter(r -> r.Modelo_Key == "d_zhou", sub_inf)[1, :]

    # Representative Sessile Models for theta0 = 40°
    sub_ses = filter(r -> r.Caso == c_id && r.theta0_deg == 40.0, df_ses_sa)
    s_ca = filter(r -> r.Modelo == "AC (CA)", sub_ses)[1, :]
    s_cr = filter(r -> r.Modelo == "RC (CR)", sub_ses)[1, :]
    s_ss = filter(r -> r.Modelo == "Stick-Slide (S-S)", sub_ses)[1, :]
    s_co = filter(r -> r.Modelo == "RC Substrate Cooled", sub_ses)[1, :]

    push!(records_summary, (
        Caso = c_id,
        T_C = t_c,
        HR_pct = hr,
        # Initial Areas [mm²]
        S0_Vuelo_mm2 = m_d.S0_mm2,
        S0_Sesil_Cap_mm2 = s_cr.S_cap0_mm2,
        # In-Flight D² Classic
        t_Vuelo_D2_s = m_a.t_evap_s,
        Smean_Vuelo_D2_mm2 = m_a.S_mean_mm2,
        Smed_Vuelo_D2_mm2  = m_a.S_median_mm2,
        # In-Flight Abramzon-Sirignano
        t_Vuelo_AS_s = m_d.t_evap_s,
        Smean_Vuelo_AS_mm2 = m_d.S_mean_mm2,
        Smed_Vuelo_AS_mm2  = m_d.S_median_mm2,
        # Sessile Mode AC (Receding)
        t_Sesil_AC_s = s_ca.t_evap_s,
        Smean_Sesil_AC_mm2 = s_ca.S_mean_mm2,
        Smed_Sesil_AC_mm2  = s_ca.S_median_mm2,
        # Sessile Mode RC (Pinned)
        t_Sesil_RC_s = s_cr.t_evap_s,
        Smean_Sesil_RC_mm2 = s_cr.S_mean_mm2,
        Smed_Sesil_RC_mm2  = s_cr.S_median_mm2,
        # Ratio of Mean Areas: Sessile RC vs In-Flight AS
        Ratio_Smean_SesilRC_vs_VueloAS = round(s_cr.S_mean_mm2 / m_d.S_mean_mm2, digits=2)
    ))
end

df_summary = DataFrame(records_summary)
summary_csv = joinpath(DATA_DIR, "table_mean_median_surface_areas_summary.csv")
CSV.write(summary_csv, df_summary)
println("Saved Master Surface Area Summary: $summary_csv")

println("\n================================================================================")
println("SURFACE AREA CALCULATIONS SUCCESSFULLY COMPLETED!")
println("================================================================================")
