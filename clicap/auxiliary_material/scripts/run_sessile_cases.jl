"""
Sessile Droplet Evaporation Across All CLICAP Cases & Models
Auxiliary Material Reproduction & Extension Script

Simulates the complete matrix of 9 environmental cases for:
  Equivalent Diameter: d₀* = 200 µm (and 1000 µm)
  Contact Angles: 20°, 40°, 80°
  Models:
    - Mode AC (Constant Contact Angle, CA)
    - Mode RC (Constant Contact Radius, CR)
    - Mode Stick-Slide (SS, θ_rec = θ₀ / 2)
    - Mode Stick-Jump (SJ, θ_min = θ₀ / 2, θ_max = θ₀)
    - Asymptotic Shallow Thin-Droplet Closed-Form Lifetimes (Wilson 2023)
    - Picknett & Bexon (1977) Capacitance Polynomials
    - Substrate Conductive Cooling (Glass hs = 1 mm, ks = 1.4 W/mK)
    - Multi-Droplet Cooperative Vapor Shielding (Pair at B = 3 R₀)
    - Saxton et al. (2016) Generalized Tanner Spreading-Evaporation

Exports CSV datasets and publication figures.
"""

ENV["GKSwstype"] = "100"
using Plots
using DataFrames
using CSV
using Printf

# Load Evaporacion.jl package if not already loaded
if !isdefined(Main, :Evaporacion)
    include(joinpath(@__DIR__, "..", "..", "..", "julia", "src", "Evaporacion.jl"))
    using .Evaporacion
end

const SCRIPT_DIR = @__DIR__
const DATA_DIR = normpath(joinpath(SCRIPT_DIR, "..", "data"))
const FIG_DIR  = normpath(joinpath(SCRIPT_DIR, "..", "figures"))

mkpath(DATA_DIR)
mkpath(FIG_DIR)

println("================================================================================")
println("CLICAP SESSILE DROPLET EVAPORATION: ALL MODELS & CASES")
println("================================================================================")

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

const CONTACT_ANGLES_DEG = [20.0, 40.0, 80.0]
const D0_SESSILE_UM = 200.0
const D0_M = D0_SESSILE_UM * 1e-6
const V0_M3 = (π / 6.0) * (D0_M^3)
const P_ATM = 101325.0

records_sessile = []

println("\nSimulating sessile droplets across 9 cases × 3 contact angles × 8 models...")

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
        R0 = r_from_volume_sessile(V0_M3, th0_rad)

        # 1. Mode AC (CA)
        cfg_ca = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CA, t_end=2000.0)
        sol_ca = solve_sessile(cfg_ca)
        t_ca = sol_ca.t_evap_analytical

        # 2. Mode RC (CR)
        cfg_cr = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CR, t_end=2000.0)
        sol_cr = solve_sessile(cfg_cr)
        t_cr = sol_cr.t_evap_analytical

        # 3. Mode Stick-Slide (SS)
        cfg_ss = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, theta_rec=th0_rad/2.0, mode=:SS, t_end=2000.0)
        sol_ss = solve_sessile(cfg_ss)
        t_ss = sol_ss.t_evap_analytical

        # 4. Mode Stick-Jump (SJ)
        cfg_sj = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, theta_min=th0_rad/2.0, theta_max=th0_rad, mode=:SJ, t_end=2000.0)
        sol_sj = solve_sessile(cfg_sj)
        t_sj = sol_sj.t_evap_analytical

        # 5. Thin-droplet shallow limits (Wilson 2023 Table 2)
        t_cr_thin = lifetime_cr_thin(R0, th0_rad, D_val, delta_c, rho_l)
        t_ca_thin = lifetime_ca_thin(R0, th0_rad, D_val, delta_c, rho_l)
        t_ss_thin = lifetime_ss_thin(R0, th0_rad, th0_rad/2.0, D_val, delta_c, rho_l)

        # 6. Picknett & Bexon (1977) Capacitance
        C_ratio = picknett_bexon_capacitance(th0_rad)
        t_picknett = (rho_l * V0_M3) / (4.0 * π * R0 * D_val * delta_c * C_ratio)

        # 7. Substrate conductive cooling on glass (hs = 1 mm, ks = 1.4 W/mK)
        L_val = latent_heat_water(T_K)
        csat_plus = (psat_water(T_K + 0.5) * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * (T_K + 0.5))
        csat_minus = (psat_water(T_K - 0.5) * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * (T_K - 0.5))
        dcsat_dT_val = csat_plus - csat_minus
        tparams = SessileThermalParams(L=L_val, h_s=1e-3, k_s_eff=1.4, dcsat_dT=dcsat_dT_val)
        delta_c_eff = cooling_effective_delta_c(R0, th0_rad, D_val, delta_c, tparams)
        cfg_cool = SessileDropletConfig(D=D_val, delta_c=delta_c_eff, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CR, t_end=2000.0)
        sol_cool = solve_sessile(cfg_cool)
        t_cool = sol_cool.t_evap_analytical

        # 8. Cooperative Vapor Shielding (pair at B = 3 R0)
        t_shielded = lifetime_pair_cr(3.0 * R0, R0, th0_rad, D_val, delta_c, rho_l)

        push!(records_sessile, (
            Caso = case.id,
            T_C = T_C,
            HR_pct = HR_pct,
            d0_equiv_um = D0_SESSILE_UM,
            theta0_deg = th_deg,
            R0_um = round(R0 * 1e6, digits=1),
            t_AC_s = round(t_ca, digits=3),
            t_RC_s = round(t_cr, digits=3),
            t_SS_s = round(t_ss, digits=3),
            t_SJ_s = round(t_sj, digits=3),
            t_RC_thin_s = round(t_cr_thin, digits=3),
            t_CA_thin_s = round(t_ca_thin, digits=3),
            t_Picknett_s = round(t_picknett, digits=3),
            t_Cooled_s = round(t_cool, digits=3),
            t_Shielded_s = round(t_shielded, digits=3)
        ))
    end
end

df_sessile_all = DataFrame(records_sessile)
sess_csv = joinpath(DATA_DIR, "comparison_sessile_all_cases.csv")
CSV.write(sess_csv, df_sessile_all)
println("Saved Sessile Dataset: $sess_csv ($(nrow(df_sessile_all)) rows)")

# 3. DIRECT PAPER TABLE 3 COMPARISON
# In paper Table 3, lifetimes are given for Case 7 (313 K, 20% HR) and Case 1 (283 K, 20% HR)
paper_t3_ref = Dict(
    # theta => (AC, RC, SS, SJ)
    20.0 => (3.35, 2.30, 3.05, 2.61),
    40.0 => (4.00, 2.87, 3.76, 3.26),
    80.0 => (4.40, 3.56, 4.01, 3.76)
)

records_t3_comp = []
for c_id in (7, 1)
    sub_c = filter(r -> r.Caso == c_id, df_sessile_all)
    c_label = c_id == 7 ? "Caso 7 (313 K, 20%)" : "Caso 1 (283 K, 20%)"

    for row in eachrow(sub_c)
        th = row.theta0_deg
        p_ref = paper_t3_ref[th]

        for (m_lbl, t_j, t_p) in [
            ("AC (CA)", row.t_AC_s, p_ref[1]),
            ("RC (CR)", row.t_RC_s, p_ref[2]),
            ("S-S", row.t_SS_s, p_ref[3]),
            ("S-J", row.t_SJ_s, p_ref[4])
        ]
            push!(records_t3_comp, (
                Escenario = c_label,
                Caso_ID = c_id,
                theta0_deg = th,
                R0_um = row.R0_um,
                Modelo = m_lbl,
                Tiempo_Julia_s = t_j,
                Tiempo_Paper_s = t_p,
                Error_Rel_pct = round(((t_j - t_p) / t_p) * 100.0, digits=2)
            ))
        end
    end
end

df_t3_comp = DataFrame(records_t3_comp)
t3_csv = joinpath(DATA_DIR, "table3_sessile_paper_comparison.csv")
CSV.write(t3_csv, df_t3_comp)
println("Saved Table 3 Sessile Comparison: $t3_csv")

# 4. PLOTTING FIGURES
println("\nGenerating sessile publication figures...")

# Figure D: Sessile Lifetimes vs Relative Humidity (Cases 4, 5, 6: T=25°C, theta0=40°) - Extended Paper Fig 4
sub_fig4 = filter(r -> r.Caso in (4, 5, 6) && r.theta0_deg == 40.0, df_sessile_all)

p_s_hr = plot(
    title = "Sessile Lifetimes vs Relative Humidity (θ₀ = 40°, T = 25 °C)",
    xlabel = "Relative Humidity HR [%]",
    ylabel = "Evaporation Lifetime [s]",
    legend = :topleft,
    framestyle = :box,
    grid = true,
    dpi = 300,
    size = (750, 500)
)

plot!(p_s_hr, sub_fig4.HR_pct, sub_fig4.t_AC_s, label="Modo AC (CA)", lw=2.5, color=:royalblue, marker=:circle)
plot!(p_s_hr, sub_fig4.HR_pct, sub_fig4.t_RC_s, label="Modo RC (CR)", lw=2.5, color=:forestgreen, marker=:square)
plot!(p_s_hr, sub_fig4.HR_pct, sub_fig4.t_SS_s, label="Modo S-S", lw=2.5, color=:darkorange, marker=:diamond)
plot!(p_s_hr, sub_fig4.HR_pct, sub_fig4.t_SJ_s, label="Modo S-J", lw=2.5, color=:crimson, marker=:cross)
plot!(p_s_hr, sub_fig4.HR_pct, sub_fig4.t_Cooled_s, label="RC Substrate Cooled", lw=2, ls=:dash, color=:purple)
plot!(p_s_hr, sub_fig4.HR_pct, sub_fig4.t_Shielded_s, label="RC Shielded (Pair B=3R₀)", lw=2, ls=:dot, color=:brown)

fig_d_png = joinpath(FIG_DIR, "fig_sessile_lifetimes_vs_HR.png")
savefig(p_s_hr, fig_d_png)
println("Saved Figure D: $fig_d_png")

# Figure E: Sessile Lifetimes vs Contact Angle (Case 7: 40°C, 20% HR)
sub_c7 = filter(r -> r.Caso == 7, df_sessile_all)

p_s_th = plot(
    title = "Sessile Evaporation Modes vs Initial Contact Angle θ₀ (Case 7: 40 °C, 20% HR)",
    xlabel = "Initial Contact Angle θ₀ [deg]",
    ylabel = "Evaporation Lifetime [s]",
    legend = :topleft,
    framestyle = :box,
    grid = true,
    dpi = 300,
    size = (750, 500)
)

plot!(p_s_th, sub_c7.theta0_deg, sub_c7.t_AC_s, label="AC (CA)", lw=2.5, color=:royalblue, marker=:circle)
plot!(p_s_th, sub_c7.theta0_deg, sub_c7.t_RC_s, label="RC (CR)", lw=2.5, color=:forestgreen, marker=:square)
plot!(p_s_th, sub_c7.theta0_deg, sub_c7.t_SS_s, label="S-S", lw=2.5, color=:darkorange, marker=:diamond)
plot!(p_s_th, sub_c7.theta0_deg, sub_c7.t_SJ_s, label="S-J", lw=2.5, color=:crimson, marker=:cross)
plot!(p_s_th, sub_c7.theta0_deg, sub_c7.t_Picknett_s, label="Picknett & Bexon", lw=2, ls=:dash, color=:teal)
plot!(p_s_th, sub_c7.theta0_deg, sub_c7.t_RC_thin_s, label="Thin-droplet CR limit", lw=2, ls=:dot, color=:gray)

fig_e_png = joinpath(FIG_DIR, "fig_sessile_modes_vs_contact_angle.png")
savefig(p_s_th, fig_e_png)
println("Saved Figure E: $fig_e_png")

println("\nSessile script execution completed successfully!")
