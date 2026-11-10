"""
Comprehensive Volume, Surface Area, V/S Ratio, and Computational Performance Benchmark
Auxiliary Material Script

1. Calculates time-averaged mean volume (V_mean), temporal median volume (V_median),
   and dynamic V/S relationships:
     V_mean / S_mean
     V_median / S_median
     Time-averaged (V/S)(t)
   for both in-flight and sessile droplets across all CLICAP environmental cases.

2. Conducts execution time benchmarking:
   - 1 initial warm-up execution per model/case (discarded to eliminate JIT compilation overhead)
   - 5 independent timed runs per model/case (high-resolution nanosecond clock)
   - Measures mean, median, min, max, and standard deviation of execution times.
   - Computes relative speedup factors across all models.

Exports publication CSV tables to clicap/auxiliary_material/data/.
"""

ENV["GKSwstype"] = "100"
using DataFrames
using CSV
using Printf
using Statistics

# Load Evaporacion.jl package if not already loaded
if !isdefined(Main, :Evaporacion)
    include(joinpath(@__DIR__, "..", "..", "..", "julia", "src", "Evaporacion.jl"))
    using .Evaporacion
end

const SCRIPT_DIR = @__DIR__
const DATA_DIR = normpath(joinpath(SCRIPT_DIR, "..", "data"))
mkpath(DATA_DIR)

println("================================================================================")
println("STARTING VOLUME, SURFACE AREA, V/S RATIO & COMPUTATIONAL TIME BENCHMARK")
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
# HELPER: Compute Mean and Median from Discrete Time Series (t, val)
# -----------------------------------------------------------------------------
function calc_mean_and_median_series(t_vec::AbstractVector{<:Real}, val_vec::AbstractVector{<:Real})
    N = length(t_vec)
    if N < 2
        return (val_vec[1], val_vec[1])
    end
    t_start = t_vec[1]
    t_end = t_vec[end]
    dt_total = t_end - t_start
    if dt_total <= 0.0
        return (val_vec[1], val_vec[1])
    end

    # Trapezoidal integration for time-averaged mean
    integral = 0.0
    for i in 1:(N - 1)
        dt = t_vec[i+1] - t_vec[i]
        integral += 0.5 * (val_vec[i] + val_vec[i+1]) * dt
    end
    val_mean = integral / dt_total

    # Median over time = value at t = t_start + dt_total / 2
    t_half = t_start + 0.5 * dt_total
    idx = searchsortedlast(t_vec, t_half)
    idx = clamp(idx, 1, N - 1)
    frac = (t_half - t_vec[idx]) / max(t_vec[idx+1] - t_vec[idx], 1e-12)
    val_median = val_vec[idx] + frac * (val_vec[idx+1] - val_vec[idx])

    return (val_mean, val_median)
end

# -----------------------------------------------------------------------------
# HELPER: Benchmark Execution Time (1 Warm-up, 5 Timed Runs)
# -----------------------------------------------------------------------------
function benchmark_fn(f::Function; n_warmup::Int=1, n_runs::Int=5)
    # 1. Warm-up run (discarded)
    for _ in 1:n_warmup
        f()
    end
    # 2. Timed runs
    times_ms = Float64[]
    for _ in 1:n_runs
        t0 = time_ns()
        f()
        dt_ms = (time_ns() - t0) * 1e-6
        push!(times_ms, dt_ms)
    end
    sort!(times_ms)
    mean_val = mean(times_ms)
    median_val = median(times_ms)
    min_val = minimum(times_ms)
    max_val = maximum(times_ms)
    std_val = n_runs > 1 ? std(times_ms) : 0.0
    return (
        mean_ms = mean_val,
        median_ms = median_val,
        min_ms = min_val,
        max_ms = max_val,
        std_ms = std_val
    )
end

# =============================================================================
# PART 1: IN-FLIGHT DROPLETS VOLUME, SURFACE AREA & V/S RATIOS
# =============================================================================
println("\n>>> [1/4] COMPUTING VOLUMES, SURFACE AREAS & V/S FOR IN-FLIGHT DROPLETS...")

records_inf = []

const INFLIGHT_MODELS = [
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
        d_mm = d_um * 1e-3
        V0_mm3 = (π / 6.0) * (d_mm^3)
        S0_mm2 = π * (d_mm^2)
        VS0_um = (V0_mm3 / S0_mm2) * 1e3 # = d_um / 6

        for (m_key, m_id, m_prof, m_desc) in INFLIGHT_MODELS
            t_evap = 0.0
            V_mean = 0.0
            V_med  = 0.0
            S_mean = 0.0
            S_med  = 0.0
            VS_mean = 0.0

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

                # Dynamic trajectories V(t), S(t), (V/S)(t)
                V_t = Float64[]
                S_t = Float64[]
                VS_t = Float64[]

                for state in sol.u
                    m_curr = max(state[1], 1e-18)
                    T_curr = state[2]
                    rho_curr = liquid_density_water(T_curr + 273.15)
                    v_curr_mm3 = (m_curr / rho_curr) * 1e9 # m^3 -> mm^3
                    d_curr_mm = cbrt(v_curr_mm3 / (π / 6.0))
                    s_curr_mm2 = π * (d_curr_mm^2)
                    vs_curr_um = (v_curr_mm3 / s_curr_mm2) * 1e3 # mm -> µm

                    push!(V_t, v_curr_mm3)
                    push!(S_t, s_curr_mm2)
                    push!(VS_t, vs_curr_um)
                end

                V_mean, V_med = calc_mean_and_median_series(sol.t, V_t)
                S_mean, S_med = calc_mean_and_median_series(sol.t, S_t)
                VS_mean, _     = calc_mean_and_median_series(sol.t, VS_t)

            elseif m_id == :vlasov_kinetic
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
                R_t_mm = [vlasov_radius_trajectory(t, R0, D_val, delta_c, rho_l0, M_w, alpha_evap, v_th) * 1e3 for t in t_grid]
                V_t = [(4.0 / 3.0) * π * (r^3) for r in R_t_mm]
                S_t = [4.0 * π * (r^2) for r in R_t_mm]
                VS_t = [(r / 3.0) * 1e3 for r in R_t_mm]

                V_mean, V_med = calc_mean_and_median_series(collect(t_grid), V_t)
                S_mean, S_med = calc_mean_and_median_series(collect(t_grid), S_t)
                VS_mean, _     = calc_mean_and_median_series(collect(t_grid), VS_t)

            elseif m_id in (:surrogate_6b, :surrogate_8dyn)
                stype = m_id == :surrogate_6b ? :modelo_6b_refit : :modelo_8_dynamic
                v_term_est = calc_terminal_velocity(d_m, T_K, P_ATM)
                cfg = SurrogateConfig(
                    model_type = stype,
                    d0_mm = d_mm,
                    T_inf = T_C,
                    T_wb = T_wb_C,
                    U = max(0.2, v_term_est),
                    HR = HR_frac,
                    T0 = T_C,
                    t_end = 2500.0
                )
                sol = solve_surrogate(cfg)
                t_evap = sol.t[end]

                V_t = [max(u[1], 1e-15) for u in sol.u]
                S_t = [π * ((6.0 * v / π)^(2.0/3.0)) for v in V_t]
                VS_t = [(v / s) * 1e3 for (v, s) in zip(V_t, S_t)]

                V_mean, V_med = calc_mean_and_median_series(sol.t, V_t)
                S_mean, S_med = calc_mean_and_median_series(sol.t, S_t)
                VS_mean, _     = calc_mean_and_median_series(sol.t, VS_t)
            end

            # Ratios
            ratio_Vmean_Smean_um = (V_mean / max(S_mean, 1e-12)) * 1e3
            ratio_Vmed_Smed_um   = (V_med  / max(S_med,  1e-12)) * 1e3

            push!(records_inf, (
                Dominio = "En Vuelo",
                Caso = case.id,
                T_C = T_C,
                HR_pct = HR_pct,
                d0_um = d_um,
                Modelo_Key = m_key,
                Modelo_Desc = m_desc,
                t_evap_s = round(t_evap, digits=4),
                V0_mm3 = round(V0_mm3, digits=8),
                S0_mm2 = round(S0_mm2, digits=6),
                VS0_um = round(VS0_um, digits=3),
                V_mean_mm3 = round(V_mean, digits=8),
                V_median_mm3 = round(V_med, digits=8),
                Ratio_Vmean_V0 = round(V_mean / V0_mm3, digits=4),
                Ratio_Vmed_V0 = round(V_med / V0_mm3, digits=4),
                S_mean_mm2 = round(S_mean, digits=6),
                S_median_mm2 = round(S_med, digits=6),
                Ratio_Smean_S0 = round(S_mean / S0_mm2, digits=4),
                Ratio_Smed_S0 = round(S_med / S0_mm2, digits=4),
                Vmean_over_Smean_um = round(ratio_Vmean_Smean_um, digits=3),
                Vmed_over_Smed_um = round(ratio_Vmed_Smed_um, digits=3),
                TimeAvg_V_over_S_um = round(VS_mean, digits=3)
            ))
        end
    end
end

df_inf_vs = DataFrame(records_inf)
CSV.write(joinpath(DATA_DIR, "volumes_surface_ratios_inflight_all_cases.csv"), df_inf_vs)
println("Saved in-flight volume & V/S dataset: $(nrow(df_inf_vs)) rows")

# =============================================================================
# PART 2: SESSILE DROPLETS VOLUME, SURFACE AREA & V/S RATIOS
# =============================================================================
println("\n>>> [2/4] COMPUTING VOLUMES, SURFACE AREAS & V/S FOR SESSILE DROPLETS...")

const CONTACT_ANGLES_DEG = [20.0, 40.0, 80.0]
const D0_SESSILE_UM = 200.0
const D0_SESSILE_M = D0_SESSILE_UM * 1e-6
const V0_M3_SES = (π / 6.0) * (D0_SESSILE_M^3)
const V0_MM3_SES = V0_M3_SES * 1e9

records_ses = []

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

        # Initial spherical cap surface area S_cap0
        S_cap0_mm2 = (2.0 * π * (R0 * 1e3)^2 / (1.0 + cos(th0_rad)))
        VS0_ses_um = (V0_MM3_SES / S_cap0_mm2) * 1e3

        # 1. Mode AC (CA)
        cfg_ca = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CA, t_end=2000.0)
        sol_ca = solve_sessile(cfg_ca)
        t_ca = sol_ca.t_evap_analytical

        # For pure AC: V(t) = V0 * (1 - t/t_evap)^(3/2), Scap(t) = Scap0 * (1 - t/t_evap)
        # Analytical integrals:
        V_mean_ca = 0.40 * V0_MM3_SES
        V_med_ca  = (1.0 / (2.0 * sqrt(2.0))) * V0_MM3_SES
        S_mean_ca = 0.50 * S_cap0_mm2
        S_med_ca  = 0.50 * S_cap0_mm2
        VS_mean_ca = (2.0 / 3.0) * VS0_ses_um

        # 2. Mode RC (CR)
        cfg_cr = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CR, t_end=2000.0)
        sol_cr = solve_sessile(cfg_cr)
        t_cr = sol_cr.t_evap_analytical

        V_cr_t = [sessile_volume_from_r(R0, th) * 1e9 for th in sol_cr.theta]
        S_cr_t = [2.0 * π * (R0 * 1e3)^2 / (1.0 + cos(th)) for th in sol_cr.theta]
        VS_cr_t = [(v / s) * 1e3 for (v, s) in zip(V_cr_t, S_cr_t)]

        V_mean_cr, V_med_cr = calc_mean_and_median_series(sol_cr.t, V_cr_t)
        S_mean_cr, S_med_cr = calc_mean_and_median_series(sol_cr.t, S_cr_t)
        VS_mean_cr, _        = calc_mean_and_median_series(sol_cr.t, VS_cr_t)

        # 3. Mode Stick-Slide (SS)
        cfg_ss = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, theta_rec=th0_rad/2.0, mode=:SS, t_end=2000.0)
        sol_ss = solve_sessile(cfg_ss)
        t_ss = sol_ss.t_evap_analytical
        V_ss_t = [sol_ss.V[i] * 1e9 for i in 1:length(sol_ss.t)]
        S_ss_t = [2.0 * π * (sol_ss.R[i] * 1e3)^2 / (1.0 + cos(sol_ss.theta[i])) for i in 1:length(sol_ss.t)]
        VS_ss_t = [(v / s) * 1e3 for (v, s) in zip(V_ss_t, S_ss_t)]

        V_mean_ss, V_med_ss = calc_mean_and_median_series(sol_ss.t, V_ss_t)
        S_mean_ss, S_med_ss = calc_mean_and_median_series(sol_ss.t, S_ss_t)
        VS_mean_ss, _        = calc_mean_and_median_series(sol_ss.t, VS_ss_t)

        # 4. Mode Stick-Jump (SJ)
        cfg_sj = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, theta_min=th0_rad/2.0, theta_max=th0_rad, mode=:SJ, t_end=2000.0)
        sol_sj = solve_sessile(cfg_sj)
        t_sj = sol_sj.t_evap_analytical
        V_sj_t = [sol_sj.V[i] * 1e9 for i in 1:length(sol_sj.t)]
        S_sj_t = [2.0 * π * (sol_sj.R[i] * 1e3)^2 / (1.0 + cos(sol_sj.theta[i])) for i in 1:length(sol_sj.t)]
        VS_sj_t = [(v / s) * 1e3 for (v, s) in zip(V_sj_t, S_sj_t)]

        V_mean_sj, V_med_sj = calc_mean_and_median_series(sol_sj.t, V_sj_t)
        S_mean_sj, S_med_sj = calc_mean_and_median_series(sol_sj.t, S_sj_t)
        VS_mean_sj, _        = calc_mean_and_median_series(sol_sj.t, VS_sj_t)

        # 5. RC Substrate Cooled
        L_val = latent_heat_water(T_K)
        csat_plus = (psat_water(T_K + 0.5) * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * (T_K + 0.5))
        csat_minus = (psat_water(T_K - 0.5) * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * (T_K - 0.5))
        dcsat_dT_val = csat_plus - csat_minus
        tparams = SessileThermalParams(L=L_val, h_s=1e-3, k_s_eff=1.4, dcsat_dT=dcsat_dT_val)
        delta_c_eff = cooling_effective_delta_c(R0, th0_rad, D_val, delta_c, tparams)
        cfg_cool = SessileDropletConfig(D=D_val, delta_c=delta_c_eff, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CR, t_end=2000.0)
        sol_cool = solve_sessile(cfg_cool)
        t_cool = sol_cool.t_evap_analytical

        V_cool_t = [sessile_volume_from_r(R0, th) * 1e9 for th in sol_cool.theta]
        S_cool_t = [2.0 * π * (R0 * 1e3)^2 / (1.0 + cos(th)) for th in sol_cool.theta]
        VS_cool_t = [(v / s) * 1e3 for (v, s) in zip(V_cool_t, S_cool_t)]

        V_mean_cool, V_med_cool = calc_mean_and_median_series(sol_cool.t, V_cool_t)
        S_mean_cool, S_med_cool = calc_mean_and_median_series(sol_cool.t, S_cool_t)
        VS_mean_cool, _          = calc_mean_and_median_series(sol_cool.t, VS_cool_t)

        # 6. RC Shielded
        t_shielded = lifetime_pair_cr(3.0 * R0, R0, th0_rad, D_val, delta_c, rho_l)
        V_mean_shield = V_mean_cr
        V_med_shield  = V_med_cr
        S_mean_shield = S_mean_cr
        S_med_shield  = S_med_cr
        VS_mean_shield = VS_mean_cr

        for (m_lbl, t_val, v_mean, v_med, s_mean, s_med, vs_mean) in [
            ("AC (CA)", t_ca, V_mean_ca, V_med_ca, S_mean_ca, S_med_ca, VS_mean_ca),
            ("RC (CR)", t_cr, V_mean_cr, V_med_cr, S_mean_cr, S_med_cr, VS_mean_cr),
            ("Stick-Slide (S-S)", t_ss, V_mean_ss, V_med_ss, S_mean_ss, S_med_ss, VS_mean_ss),
            ("Stick-Jump (S-J)", t_sj, V_mean_sj, V_med_sj, S_mean_sj, S_med_sj, VS_mean_sj),
            ("RC Substrate Cooled", t_cool, V_mean_cool, V_med_cool, S_mean_cool, S_med_cool, VS_mean_cool),
            ("RC Shielded (B=3R₀)", t_shielded, V_mean_shield, V_med_shield, S_mean_shield, S_med_shield, VS_mean_shield)
        ]
            push!(records_ses, (
                Dominio = "Sésil",
                Caso = case.id,
                T_C = T_C,
                HR_pct = HR_pct,
                d0_equiv_um = D0_SESSILE_UM,
                theta0_deg = th_deg,
                R0_um = round(R0 * 1e6, digits=1),
                Modelo = m_lbl,
                t_evap_s = round(t_val, digits=3),
                V0_mm3 = round(V0_MM3_SES, digits=8),
                S_cap0_mm2 = round(S_cap0_mm2, digits=6),
                VS0_um = round(VS0_ses_um, digits=3),
                V_mean_mm3 = round(v_mean, digits=8),
                V_median_mm3 = round(v_med, digits=8),
                Ratio_Vmean_V0 = round(v_mean / V0_MM3_SES, digits=4),
                Ratio_Vmed_V0 = round(v_med / V0_MM3_SES, digits=4),
                S_mean_mm2 = round(s_mean, digits=6),
                S_median_mm2 = round(s_med, digits=6),
                Ratio_Smean_S0 = round(s_mean / S_cap0_mm2, digits=4),
                Ratio_Smed_S0 = round(s_med / S_cap0_mm2, digits=4),
                Vmean_over_Smean_um = round((v_mean / s_mean) * 1e3, digits=3),
                Vmed_over_Smed_um = round((v_med / s_med) * 1e3, digits=3),
                TimeAvg_V_over_S_um = round(vs_mean, digits=3)
            ))
        end
    end
end

df_ses_vs = DataFrame(records_ses)
CSV.write(joinpath(DATA_DIR, "volumes_surface_ratios_sessile_all_cases.csv"), df_ses_vs)
println("Saved sessile volume & V/S dataset: $(nrow(df_ses_vs)) rows")

# =============================================================================
# PART 3: MASTER CROSS-DOMAIN VOLUME & V/S SUMMARY (d0 = 200 µm, θ0 = 40°)
# =============================================================================
println("\n>>> [3/4] GENERATING MASTER V/S COMPARATIVE SUMMARY TABLE...")

records_master_vs = []

for c_id in 1:9
    sub_inf = filter(r -> r.Caso == c_id && r.d0_um == 200.0, df_inf_vs)
    sub_ses = filter(r -> r.Caso == c_id && r.theta0_deg == 40.0, df_ses_vs)

    r_as = filter(r -> r.Modelo_Key == "d", sub_inf)[1, :]
    r_d2 = filter(r -> r.Modelo_Key == "a", sub_inf)[1, :]
    r_ca = filter(r -> r.Modelo == "AC (CA)", sub_ses)[1, :]
    r_cr = filter(r -> r.Modelo == "RC (CR)", sub_ses)[1, :]

    push!(records_master_vs, (
        Caso = c_id,
        T_C = r_as.T_C,
        HR_pct = r_as.HR_pct,
        V0_mm3 = r_as.V0_mm3,
        S0_Vuelo_mm2 = r_as.S0_mm2,
        S0_Sesil_Cap_mm2 = r_cr.S_cap0_mm2,
        VS0_Vuelo_um = r_as.VS0_um,
        VS0_Sesil_um = r_cr.VS0_um,
        # In-Flight A-S
        Vmean_AS_mm3 = r_as.V_mean_mm3,
        Smean_AS_mm2 = r_as.S_mean_mm2,
        VSmean_AS_um = r_as.Vmean_over_Smean_um,
        # Sessile AC (Receding)
        Vmean_SesilAC_mm3 = r_ca.V_mean_mm3,
        Smean_SesilAC_mm2 = r_ca.S_mean_mm2,
        VSmean_SesilAC_um = r_ca.Vmean_over_Smean_um,
        # Sessile RC (Pinned)
        Vmean_SesilRC_mm3 = r_cr.V_mean_mm3,
        Smean_SesilRC_mm2 = r_cr.S_mean_mm2,
        VSmean_SesilRC_um = r_cr.Vmean_over_Smean_um,
        # Cross-Domain Ratios
        Ratio_VSmean_RC_vs_AS = round(r_cr.Vmean_over_Smean_um / r_as.Vmean_over_Smean_um, digits=3),
        Ratio_Smean_RC_vs_AS = round(r_cr.S_mean_mm2 / r_as.S_mean_mm2, digits=3)
    ))
end

df_master_vs = DataFrame(records_master_vs)
CSV.write(joinpath(DATA_DIR, "table_mean_median_volume_surface_summary.csv"), df_master_vs)
println("Saved Master Volume & V/S Summary: $(nrow(df_master_vs)) rows")

# =============================================================================
# PART 4: COMPUTATIONAL TIME ANALYSIS (1 Warm-up, 5 Timed Runs)
# =============================================================================
println("\n>>> [4/4] BENCHMARKING MODEL CALCULATION TIMES (1 Warmup Discarded, 5 Runs)...")

records_timings = []

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

    d_um = 200.0
    d_m = d_um * 1e-6
    d_mm = d_um * 1e-3
    rho_l0 = liquid_density_water(T_C + 273.15)
    m0 = rho_l0 * (π / 6.0) * (d_m^3)
    tspan = (0.0, 100.0)

    # 1. In-flight models
    for (m_key, m_id, m_prof, m_desc) in INFLIGHT_MODELS
        bench = if m_id in (:d2_classic, :d2_stefan, :d2_convective, :d2_revised_dallabarba, :abramzon_sirignano)
            isothermal = m_id in (:d2_classic, :d2_stefan, :d2_convective, :d2_revised_dallabarba)
            is_conv = m_id in (:d2_convective, :d2_revised_dallabarba, :abramzon_sirignano)
            T_init = isothermal ? T_wb_C : T_C
            p = DropletParams(
                T_inf = T_C, Y_inf = Y_inf, P_atm = P_ATM,
                model_variant = m_id, profile = m_prof,
                terminal_velocity = is_conv, isothermal = isothermal, T_d_fixed = T_wb_C
            )
            benchmark_fn(() -> solve_droplet(p, [m0, T_init], tspan); n_warmup=1, n_runs=5)

        elseif m_id == :vlasov_kinetic
            D_val = gas_diffusivity_water_air(T_K, P_ATM)
            delta_c = (psat_water(T_wb_C + 273.15) / (Evaporacion.Common.Constants.R_GAS * (T_wb_C + 273.15)) * M_w) * (1.0 - HR_frac)
            R0 = d_m / 2.0
            alpha_evap = 0.04
            v_th = sqrt(8.0 * Evaporacion.Common.Constants.R_GAS * T_K / (π * M_w))
            benchmark_fn(() -> begin
                tau_diff = (rho_l0 * R0^2) / (2.0 * D_val * delta_c)
                l_k = vlasov_kinetic_length(D_val, alpha_evap, v_th)
                tau_diff * (1.0 + 2.0 * l_k / R0)
            end; n_warmup=1, n_runs=5)

        elseif m_id in (:surrogate_6b, :surrogate_8dyn)
            stype = m_id == :surrogate_6b ? :modelo_6b_refit : :modelo_8_dynamic
            v_term = calc_terminal_velocity(d_m, T_K, P_ATM)
            cfg_s = SurrogateConfig(
                model_type = stype, d0_mm = d_mm, T_inf = T_C, T_wb = T_wb_C,
                U = max(0.2, v_term), HR = HR_frac, T0 = T_C, t_end = 2500.0
            )
            benchmark_fn(() -> solve_surrogate(cfg_s); n_warmup=1, n_runs=5)
        end

        push!(records_timings, (
            Dominio = "En Vuelo",
            Caso = case.id,
            T_C = T_C,
            HR_pct = HR_pct,
            Modelo_Key = m_key,
            Modelo_Desc = m_desc,
            Mean_ms = round(bench.mean_ms, digits=4),
            Median_ms = round(bench.median_ms, digits=4),
            Min_ms = round(bench.min_ms, digits=4),
            Max_ms = round(bench.max_ms, digits=4),
            Std_ms = round(bench.std_ms, digits=4)
        ))
    end

    # 2. In-flight 1D Radial Moving-Boundary PDE (N = 11)
    cfg_pde = RadialPDEConfig(
        n_nodes = 11, d0 = d_m, T_drop_0 = T_C, T_air = T_C, RH = HR_frac, t_end = 10.0
    )
    bench_pde = benchmark_fn(() -> solve_radial_pde(cfg_pde); n_warmup=1, n_runs=5)
    push!(records_timings, (
        Dominio = "En Vuelo",
        Caso = case.id,
        T_C = T_C,
        HR_pct = HR_pct,
        Modelo_Key = "pde",
        Modelo_Desc = "pde (1D Radial Conduction PDE, N=11)",
        Mean_ms = round(bench_pde.mean_ms, digits=4),
        Median_ms = round(bench_pde.median_ms, digits=4),
        Min_ms = round(bench_pde.min_ms, digits=4),
        Max_ms = round(bench_pde.max_ms, digits=4),
        Std_ms = round(bench_pde.std_ms, digits=4)
    ))

    # 3. Sessile droplet models (theta0 = 40°)
    th0_rad = deg2rad(40.0)
    R0 = r_from_volume_sessile(V0_M3_SES, th0_rad)
    D_val = gas_diffusivity_water_air(T_K, P_ATM)
    c_sat = (psat_water(T_K) * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * T_K)
    delta_c = (1.0 - HR_frac) * c_sat

    # Sessile AC (CA)
    cfg_ca = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l0, R0=R0, theta0=th0_rad, mode=:CA, t_end=2000.0)
    bench_ca = benchmark_fn(() -> solve_sessile(cfg_ca); n_warmup=1, n_runs=5)
    push!(records_timings, (
        Dominio = "Sésil", Caso = case.id, T_C = T_C, HR_pct = HR_pct,
        Modelo_Key = "AC", Modelo_Desc = "AC (CA - Constant Contact Angle)",
        Mean_ms = round(bench_ca.mean_ms, digits=4), Median_ms = round(bench_ca.median_ms, digits=4),
        Min_ms = round(bench_ca.min_ms, digits=4), Max_ms = round(bench_ca.max_ms, digits=4),
        Std_ms = round(bench_ca.std_ms, digits=4)
    ))

    # Sessile RC (CR)
    cfg_cr = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l0, R0=R0, theta0=th0_rad, mode=:CR, t_end=2000.0)
    bench_cr = benchmark_fn(() -> solve_sessile(cfg_cr); n_warmup=1, n_runs=5)
    push!(records_timings, (
        Dominio = "Sésil", Caso = case.id, T_C = T_C, HR_pct = HR_pct,
        Modelo_Key = "RC", Modelo_Desc = "RC (CR - Constant Contact Radius)",
        Mean_ms = round(bench_cr.mean_ms, digits=4), Median_ms = round(bench_cr.median_ms, digits=4),
        Min_ms = round(bench_cr.min_ms, digits=4), Max_ms = round(bench_cr.max_ms, digits=4),
        Std_ms = round(bench_cr.std_ms, digits=4)
    ))

    # Sessile Stick-Slide (SS)
    cfg_ss = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l0, R0=R0, theta0=th0_rad, theta_rec=th0_rad/2.0, mode=:SS, t_end=2000.0)
    bench_ss = benchmark_fn(() -> solve_sessile(cfg_ss); n_warmup=1, n_runs=5)
    push!(records_timings, (
        Dominio = "Sésil", Caso = case.id, T_C = T_C, HR_pct = HR_pct,
        Modelo_Key = "SS", Modelo_Desc = "SS (Stick-Slide)",
        Mean_ms = round(bench_ss.mean_ms, digits=4), Median_ms = round(bench_ss.median_ms, digits=4),
        Min_ms = round(bench_ss.min_ms, digits=4), Max_ms = round(bench_ss.max_ms, digits=4),
        Std_ms = round(bench_ss.std_ms, digits=4)
    ))

    # Sessile Stick-Jump (SJ)
    cfg_sj = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l0, R0=R0, theta0=th0_rad, theta_min=th0_rad/2.0, theta_max=th0_rad, mode=:SJ, t_end=2000.0)
    bench_sj = benchmark_fn(() -> solve_sessile(cfg_sj); n_warmup=1, n_runs=5)
    push!(records_timings, (
        Dominio = "Sésil", Caso = case.id, T_C = T_C, HR_pct = HR_pct,
        Modelo_Key = "SJ", Modelo_Desc = "SJ (Stick-Jump)",
        Mean_ms = round(bench_sj.mean_ms, digits=4), Median_ms = round(bench_sj.median_ms, digits=4),
        Min_ms = round(bench_sj.min_ms, digits=4), Max_ms = round(bench_sj.max_ms, digits=4),
        Std_ms = round(bench_sj.std_ms, digits=4)
    ))
end

df_timings_all = DataFrame(records_timings)
CSV.write(joinpath(DATA_DIR, "calculation_times_all_cases.csv"), df_timings_all)
println("Saved detailed calculation timings: $(nrow(df_timings_all)) rows")

# Compute overall average timings and speedup relative to Abramzon-Sirignano
grouped = groupby(df_timings_all, [:Dominio, :Modelo_Key, :Modelo_Desc])
summary_timings = combine(grouped,
    :Mean_ms => mean => :Mean_ms,
    :Median_ms => mean => :Median_ms,
    :Min_ms => minimum => :Min_ms,
    :Max_ms => maximum => :Max_ms,
    :Std_ms => mean => :Std_ms
)

# Reference timing = Abramzon-Sirignano (d)
t_as_ref = filter(r -> r.Modelo_Key == "d", summary_timings).Mean_ms[1]
summary_timings[!, :Speedup_vs_AS] = round.(t_as_ref ./ summary_timings.Mean_ms, digits=2)

sort!(summary_timings, :Mean_ms)
CSV.write(joinpath(DATA_DIR, "table_model_computation_times_summary.csv"), summary_timings)
println("Saved Summary Calculation Times Table: $(nrow(summary_timings)) models")

println("\n================================================================================")
println("CALCULATION TIME SUMMARY TABLE (Mean over 9 Environmental Cases, 5 Runs Each)")
println("================================================================================")
for r in eachrow(summary_timings)
    @printf("%-10s | %-12s | Mean: %8.4f ms | Med: %8.4f ms | Min: %8.4f ms | Max: %8.4f ms | Speedup vs A-S: %6.2f×\n",
        r.Dominio, r.Modelo_Key, r.Mean_ms, r.Median_ms, r.Min_ms, r.Max_ms, r.Speedup_vs_AS)
end
println("================================================================================")
