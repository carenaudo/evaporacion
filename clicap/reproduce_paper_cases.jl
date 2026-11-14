"""
Paper Reproduction Script (Renaudo, Herrera, Bertin - PLAPIQUI / UNS / CONICET)
Executes all models and scenarios from Section 2.4 strictly in Julia:
  - Table 1: Matrix of 9 Environmental Cases
  - Table 2: Lifetimes of In-Flight Droplets (Models a, b, c, d for Cases 1 and 7)
  - Figure 3: Initial Evaporation Rates K(t=0) for 200 µm droplets (Cases 4, 5, 6)
  - Table 3: Lifetimes of Sessile Droplets (AC, RC, S-S, S-J for theta0 = 20°, 40°, 80°)
  - Figure 4: Sessile Droplet Lifetimes vs Relative Humidity (Cases 4, 5, 6)
Exports corresponding CSV tables and PNG figures.
"""

using Plots
using DataFrames
using CSV
using OrdinaryDiffEq

include(joinpath(@__DIR__, "..", "julia", "src", "Evaporacion.jl"))
using .Evaporacion

const ARTIFACT_DIR = @__DIR__
mkpath(ARTIFACT_DIR)

println("================================================================================")
println("REPRODUCING PAPER: Comparación de modelos de evaporación (Renaudo et al.)")
println("================================================================================")

# ==============================================================================
# 1. TABLE 1: ENVIRONMENTAL CASES
# ==============================================================================
println("\n1. Generating Table 1 (Environmental Cases)...")

cases_data = [
    (1, 283.15, 20.0),
    (2, 283.15, 50.0),
    (3, 283.15, 80.0),
    (4, 298.15, 20.0),
    (5, 298.15, 50.0),
    (6, 298.15, 80.0),
    (7, 313.15, 20.0),
    (8, 313.15, 50.0),
    (9, 313.15, 80.0)
]

df_table1 = DataFrame(
    Caso = [c[1] for c in cases_data],
    T_K = [round(c[2], digits=0) for c in cases_data],
    T_C = [round(c[2] - 273.15, digits=1) for c in cases_data],
    HR_pct = [round(c[3], digits=0) for c in cases_data]
)

table1_csv = joinpath(ARTIFACT_DIR, "table1_environmental_cases.csv")
CSV.write(table1_csv, df_table1)
println("   Saved Table 1 CSV: $table1_csv")

# ==============================================================================
# 2. TABLE 2: IN-FLIGHT DROPLETS (CASES 1 AND 7)
# ==============================================================================
println("\n2. Simulating In-Flight Droplet Lifetimes (Table 2)...")

diameters_um = [50.0, 200.0, 600.0, 1200.0]
diameters_m = diameters_um .* 1e-6
models = [:d2_classic, :d2_stefan, :d2_convective, :abramzon_sirignano]
model_labels = ["a (Ley D² clásica)", "b (Ley D² c/ Stefan)", "c (Ley D² c/ Convección)", "d (Abramzon-Sirignano)"]

# Container for Table 2 results
records_t2 = []

# Table 2 of the RevICAp 2026 draft, for comparison. These were computed with an absolute
# tolerance of 1e-9 kg on the mass (larger than the whole mass of the smaller drops) and with
# model d using the parabolic profile; replace them with this script's output once rerun.
paper_ref = Dict(
    # (model, d0_um) => (Case 1 draft, Case 7 draft)
    (:d2_classic, 50.0) => (3.46, 1.29),
    (:d2_stefan, 50.0) => (3.44, 1.26),
    (:d2_convective, 50.0) => (3.24, 1.17),
    (:abramzon_sirignano, 50.0) => (2.9, 1.0),

    (:d2_classic, 200.0) => (78.66, 25.97),
    (:d2_stefan, 200.0) => (78.43, 25.63),
    (:d2_convective, 200.0) => (32.45, 17.68),
    (:abramzon_sirignano, 200.0) => (34.99, 12.48),

    (:d2_classic, 600.0) => (535.65, 203.78),
    (:d2_stefan, 600.0) => (532.28, 192.36),
    (:d2_convective, 600.0) => (196.41, 94.5),
    (:abramzon_sirignano, 600.0) => (189.99, 69.62),

    (:d2_classic, 1200.0) => (2143.87, 738.52),
    (:d2_stefan, 1200.0) => (2137.95, 734.49),
    (:d2_convective, 1200.0) => (616.29, 218.97),
    (:abramzon_sirignano, 1200.0) => (495.14, 183.96)
)

function simulate_inflight_case(m_sym, d0, T_K, HR_frac)
    T_inf_C = T_K - 273.15
    P_atm = 101325.0
    psat_inf = psat_water(T_K)
    pv_inf = HR_frac * psat_inf
    M_w = Evaporacion.Common.Constants.M_WATER
    M_a = Evaporacion.Common.Constants.M_AIR
    Y_inf = (pv_inf * M_w) / (pv_inf * M_w + (P_atm - pv_inf) * M_a)
    T_wb_C = stull_twb(T_inf_C, HR_frac * 100.0)

    # Initial droplet mass at wet bulb (for models a, b, c) or at T_inf (for model d)
    T_init_C = m_sym == :abramzon_sirignano ? T_inf_C : T_wb_C
    rho_l0 = liquid_density_water(T_init_C + 273.15)
    m0 = rho_l0 * (π / 6.0) * (d0^3)

    isothermal = m_sym !== :abramzon_sirignano
    is_conv = m_sym in (:d2_convective, :abramzon_sirignano)

    params = DropletParams(
        T_inf = T_inf_C,
        Y_inf = Y_inf,
        P_atm = P_atm,
        model_variant = m_sym,
        terminal_velocity = is_conv,
        isothermal = isothermal,
        T_d_fixed = T_wb_C,
        profile = :uniforme,     # model d: uniform droplet temperature (Methods, d.3)
        frossling_coef = 0.553   # Eqs. 8 and 13
    )

    t_est = (d0 / 50e-6)^2 * 10.0
    tspan = (0.0, max(100.0, t_est * 3.0))

    sol = solve_droplet(params, [m0, T_init_C], tspan)
    # sol.t[end] is the lifetime only if the mass callback stopped the run
    string(sol.retcode) == "Terminated" ||
        error("$m_sym, d0 = $(d0 * 1e6) µm, T = $T_K K did not evaporate within tspan ($(sol.retcode))")
    return sol.t[end]
end

for d_um in diameters_um
    d_m = d_um * 1e-6
    for (m_idx, m_sym) in enumerate(models)
        # Case 1 (283K, 20% HR)
        t_c1_true = simulate_inflight_case(m_sym, d_m, 283.15, 0.20)
        # Case 7 (313K, 20% HR)
        t_c7_true = simulate_inflight_case(m_sym, d_m, 313.15, 0.20)

        # Paper reported values
        p_c1, p_c7 = paper_ref[(m_sym, d_um)]

        push!(records_t2, (
            Modelo = ['a', 'b', 'c', 'd'][m_idx],
            Descripcion = model_labels[m_idx],
            d0_um = d_um,
            t_Caso1_Julia_s = round(t_c1_true, digits=2),
            t_Caso7_Julia_s = round(t_c7_true, digits=2),
            t_Caso1_Paper_s = p_c1,
            t_Caso7_Paper_s = p_c7
        ))
    end
end

df_table2 = DataFrame(records_t2)
table2_csv = joinpath(ARTIFACT_DIR, "table2_inflight_lifetimes.csv")
CSV.write(table2_csv, df_table2)
println("   Saved Table 2 CSV: $table2_csv")

# ==============================================================================
# 3. FIGURE 3: INITIAL EVAPORATION COEFFICIENT K(t=0)
# ==============================================================================
println("\n3. Computing Figure 3: Initial Evaporation Rates K(t=0) for d0 = 200 µm...")

cases_fig3 = [(4, 298.15, 0.20), (5, 298.15, 0.50), (6, 298.15, 0.80)]
d0_f3 = 200e-6

records_f3 = []

for (c_num, T_K, HR_val) in cases_fig3
    T_inf_C = T_K - 273.15
    P_atm = 101325.0
    psat_inf = psat_water(T_K)
    pv_inf = HR_val * psat_inf
    M_w = Evaporacion.Common.Constants.M_WATER
    M_a = Evaporacion.Common.Constants.M_AIR
    Y_inf = (pv_inf * M_w) / (pv_inf * M_w + (P_atm - pv_inf) * M_a)
    T_wb_C = stull_twb(T_inf_C, HR_val * 100.0)

    for (m_idx, m_sym) in enumerate(models)
        T_init_C = m_sym == :abramzon_sirignano ? T_inf_C : T_wb_C
        rho_l0 = liquid_density_water(T_init_C + 273.15)
        m0 = rho_l0 * (π / 6.0) * (d0_f3^3)

        isothermal = m_sym !== :abramzon_sirignano
        is_conv = m_sym in (:d2_convective, :abramzon_sirignano)

        params = DropletParams(
            T_inf = T_inf_C,
            Y_inf = Y_inf,
            P_atm = P_atm,
            model_variant = m_sym,
            terminal_velocity = is_conv,
            isothermal = isothermal,
            T_d_fixed = T_wb_C,
            profile = :uniforme,
            frossling_coef = 0.553
        )

        # Evaluate initial algebraic state at t=0
        alg0 = calc_droplet_algebraic(m0, T_init_C, T_init_C, params)
        # K = 4 * m_dot / (pi * rho_l * d) in mm^2/s
        K_m2_s = 4.0 * alg0.m_dot / (π * alg0.rho_l * alg0.d)
        K_mm2_s = K_m2_s * 1e6

        push!(records_f3, (
            Caso = c_num,
            HR_pct = Int(HR_val * 100),
            Modelo = ['a', 'b', 'c', 'd'][m_idx],
            Modelo_Nombre = ["a (Clasica)", "b (Stefan)", "c (Conveccion)", "d (A-S)"][m_idx],
            K_t0_mm2_s = round(K_mm2_s, digits=5)
        ))
    end
end

df_fig3 = DataFrame(records_f3)
fig3_csv = joinpath(ARTIFACT_DIR, "figure3_initial_evaporation_rates.csv")
CSV.write(fig3_csv, df_fig3)
println("   Saved Figure 3 CSV: $fig3_csv")

# Plot Figure 3 as a grouped bar chart
p_fig3 = plot(
    title = "Figura 3: Coeficiente de evaporación a t=0 (d0 = 200 µm, T = 298 K)",
    xlabel = "Humedad Relativa (HR)",
    ylabel = "K(t=0) [mm²/s]",
    framestyle = :box,
    grid = true,
    fontfamily = "Times",
    titlefontsize = 14,
    guidefontsize = 14,
    tickfontsize = 12,
    legendfontsize = 11,
    left_margin = 5Plots.mm,
    bottom_margin = 4Plots.mm,
    ylim = (0.0, maximum(df_fig3.K_t0_mm2_s) * 1.12),
    dpi = 300,
    size = (750, 500)
)

categories = ["HR = 20% (Caso 4)", "HR = 50% (Caso 5)", "HR = 80% (Caso 6)"]
x_pos = [1.0, 2.0, 3.0]
bar_width = 0.18
offsets = [-1.5 * bar_width, -0.5 * bar_width, 0.5 * bar_width, 1.5 * bar_width]
colors = [:royalblue, :forestgreen, :darkorange, :crimson]

for (m_idx, m_char) in enumerate(['a', 'b', 'c', 'd'])
    sub_df = filter(row -> row.Modelo == m_char, df_fig3)
    y_vals = sub_df.K_t0_mm2_s
    bar!(
        p_fig3,
        x_pos .+ offsets[m_idx],
        y_vals,
        bar_width = bar_width,
        label = ["a (D² Clásica)", "b (D² Stefan)", "c (D² Convección)", "d (Abramzon-Sirignano)"][m_idx],
        color = colors[m_idx],
        alpha = 0.85
    )
end

xticks!(p_fig3, x_pos, categories)
fig3_png = joinpath(ARTIFACT_DIR, "figure3_initial_evaporation_rates.png")
savefig(p_fig3, fig3_png)
println("   Saved Figure 3 PNG: $fig3_png")

# ==============================================================================
# 4. TABLE 3: SESSILE DROPLETS (AC, RC, S-S, S-J)
# ==============================================================================
println("\n4. Simulating Sessile Droplet Lifetimes (Table 3)...")

d0_sessile = 200e-6
V0_sessile = sphere_volume_from_d(d0_sessile)
contact_angles_deg = [20.0, 40.0, 80.0]

# Table 3 of the RevICAp 2026 draft (Case 7), for comparison
paper_t3_ref = Dict(
    # theta => (AC, RC, SS, SJ)
    20.0 => (5.32, 3.66, 4.54, 4.86),
    40.0 => (6.34, 4.54, 5.55, 5.89),
    80.0 => (7.03, 5.65, 6.56, 6.81)
)

records_t3 = []

# We evaluate Table 3 for Case 7 (313K, 20% HR) and Case 1 (283K, 20% HR)
for T_eval in [313.15, 283.15]
    c_name = T_eval == 313.15 ? "Caso 7 (313 K)" : "Caso 1 (283 K)"
    rho_l = liquid_density_water(T_eval)
    D_val = gas_diffusivity_water_air(T_eval, 101325.0)
    psat = psat_water(T_eval)
    c_sat = (psat * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * T_eval)
    delta_c = 0.80 * c_sat  # HR = 20% -> (1 - 0.20) * csat

    for th_deg in contact_angles_deg
        th0 = deg2rad(th_deg)
        R0 = r_from_volume_sessile(V0_sessile, th0)

        # AC
        cfg_ca = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0, mode=:CA, t_end=200.0)
        sol_ca = solve_sessile(cfg_ca)

        # RC
        cfg_cr = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0, mode=:CR, t_end=200.0)
        sol_cr = solve_sessile(cfg_cr)

        # S-S (theta_rec = theta0 / 2)
        cfg_ss = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0, theta_rec=th0/2.0, mode=:SS, t_end=200.0)
        sol_ss = solve_sessile(cfg_ss)

        # S-J (theta_min = theta0 / 2, theta_max = theta0)
        cfg_sj = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0, theta_min=th0/2.0, theta_max=th0, mode=:SJ, t_end=200.0)
        sol_sj = solve_sessile(cfg_sj)

        # the draft reports Case 7 only
        p_vals = T_eval == 313.15 ? paper_t3_ref[th_deg] : (missing, missing, missing, missing)

        for (m_label, t_ana, p_val) in [
            ("AC", sol_ca.t_evap_analytical, p_vals[1]),
            ("RC", sol_cr.t_evap_analytical, p_vals[2]),
            ("S-S", sol_ss.t_evap_analytical, p_vals[3]),
            ("S-J", sol_sj.t_evap_analytical, p_vals[4])
        ]
            push!(records_t3, (
                Escenario = c_name,
                Modelo = m_label,
                d0_star_um = 200.0,
                AC_deg = th_deg,
                R0_um = round(R0 * 1e6, digits=1),
                Tiempo_vida_Julia_s = round(t_ana, digits=2),
                Tiempo_vida_Paper_s = p_val
            ))
        end
    end
end

df_table3 = DataFrame(records_t3)
table3_csv = joinpath(ARTIFACT_DIR, "table3_sessile_lifetimes.csv")
CSV.write(table3_csv, df_table3)
println("   Saved Table 3 CSV: $table3_csv")

# ==============================================================================
# 5. FIGURE 4: SESSILE DROPLET LIFETIMES VS RELATIVE HUMIDITY (CASES 4, 5, 6)
# ==============================================================================
println("\n5. Generating Figure 4: Sessile Droplet Lifetimes for Cases 4, 5, 6 (T = 298 K)...")

records_f4 = []
cases_f4 = [(4, 20.0), (5, 50.0), (6, 80.0)]
T_f4 = 298.15
rho_l_f4 = liquid_density_water(T_f4)
D_f4 = gas_diffusivity_water_air(T_f4, 101325.0)
psat_f4 = psat_water(T_f4)
c_sat_f4 = (psat_f4 * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * T_f4)

for (c_num, hr_pct) in cases_f4
    delta_c_f4 = (1.0 - hr_pct / 100.0) * c_sat_f4
    for th_deg in [20.0, 40.0, 80.0]
        th0 = deg2rad(th_deg)
        R0 = r_from_volume_sessile(V0_sessile, th0)

        cfg_ca = SessileDropletConfig(D=D_f4, delta_c=delta_c_f4, rho_l=rho_l_f4, R0=R0, theta0=th0, mode=:CA, t_end=500.0)
        sol_ca = solve_sessile(cfg_ca)

        cfg_cr = SessileDropletConfig(D=D_f4, delta_c=delta_c_f4, rho_l=rho_l_f4, R0=R0, theta0=th0, mode=:CR, t_end=500.0)
        sol_cr = solve_sessile(cfg_cr)

        cfg_ss = SessileDropletConfig(D=D_f4, delta_c=delta_c_f4, rho_l=rho_l_f4, R0=R0, theta0=th0, theta_rec=th0/2.0, mode=:SS, t_end=500.0)
        sol_ss = solve_sessile(cfg_ss)

        cfg_sj = SessileDropletConfig(D=D_f4, delta_c=delta_c_f4, rho_l=rho_l_f4, R0=R0, theta0=th0, theta_min=th0/2.0, theta_max=th0, mode=:SJ, t_end=500.0)
        sol_sj = solve_sessile(cfg_sj)

        push!(records_f4, (Caso=c_num, HR_pct=hr_pct, theta0_deg=th_deg, Modo="AC", t_evap_s=round(sol_ca.t_evap_analytical, digits=2)))
        push!(records_f4, (Caso=c_num, HR_pct=hr_pct, theta0_deg=th_deg, Modo="RC", t_evap_s=round(sol_cr.t_evap_analytical, digits=2)))
        push!(records_f4, (Caso=c_num, HR_pct=hr_pct, theta0_deg=th_deg, Modo="S-S", t_evap_s=round(sol_ss.t_evap_analytical, digits=2)))
        push!(records_f4, (Caso=c_num, HR_pct=hr_pct, theta0_deg=th_deg, Modo="S-J", t_evap_s=round(sol_sj.t_evap_analytical, digits=2)))
    end
end

df_fig4 = DataFrame(records_f4)
fig4_csv = joinpath(ARTIFACT_DIR, "figure4_sessile_lifetimes.csv")
CSV.write(fig4_csv, df_fig4)
println("   Saved Figure 4 CSV: $fig4_csv")

# Plot Figure 4: Lifetime vs HR for each mode (at theta0 = 40°) as a grouped bar chart
sub_theta40 = filter(r -> r.theta0_deg == 40.0, df_fig4)
p_fig4 = plot(
    title = "Figura 4: Tiempo de vida sésil (d0* = 200 µm, θ0 = 40°, T = 298 K)",
    xlabel = "Humedad Relativa (HR)",
    ylabel = "Tiempo de vida [s]",
    framestyle = :box,
    grid = true,
    fontfamily = "Times",
    titlefontsize = 14,
    guidefontsize = 14,
    tickfontsize = 12,
    legendfontsize = 11,
    left_margin = 5Plots.mm,
    bottom_margin = 4Plots.mm,
    ylim = (0.0, maximum(sub_theta40.t_evap_s) * 1.12),
    dpi = 300,
    size = (750, 500),
    legend = :topleft
)

modes_list = ["AC", "RC", "S-S", "S-J"]
mode_colors = Dict("AC" => :forestgreen, "RC" => :navy, "S-S" => :purple, "S-J" => :darkorange)

for (m_idx, m) in enumerate(modes_list)
    sub = sort(filter(r -> r.Modo == m, sub_theta40), :HR_pct)
    bar!(
        p_fig4,
        x_pos .+ offsets[m_idx],
        sub.t_evap_s,
        bar_width = bar_width,
        label = "Modo $m (θ0 = 40°)",
        color = mode_colors[m],
        alpha = 0.85
    )
end

xticks!(p_fig4, x_pos, categories)
fig4_png = joinpath(ARTIFACT_DIR, "figure4_sessile_lifetimes.png")
savefig(p_fig4, fig4_png)
println("   Saved Figure 4 PNG: $fig4_png")

println("\n================================================================================")
println("REPRODUCTION SIMULATION COMPLETED SUCCESSFULLY!")
println("================================================================================")
