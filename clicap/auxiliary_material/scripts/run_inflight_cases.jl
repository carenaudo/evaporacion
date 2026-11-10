"""
In-Flight Droplet Evaporation Across All CLICAP Cases & Models
Auxiliary Material Reproduction & Extension Script

Simulates the complete matrix of 9 environmental cases for:
  Diameters: 50, 200, 600, 1200 µm
  Models:
    - Model a: D² classic (Langmuir/Godsave, Sh = 2, Twb)
    - Model b: D² with Stefan blowing (Spalding, ln(1+BM), Twb)
    - Model c: D² with dynamic convection (Frössling, v_term(d), Twb)
    - Model d: Abramzon & Sirignano (1989) parabolic profile
    - Model d_uni: Abramzon & Sirignano uniform profile
    - Model d_zhou: Zhou et al. (2017) 3rd-order profile with Hill's vortex stagnation
    - Model revised: Dalla Barba et al. (2021) revised D² with ventilated Twb(Re)
    - Model vlasov: Vlasov (2021) kinetic-diffusive Knudsen transition
    - Model surr_6b: Calibrated surrogate Modelo 6b Refit
    - Model surr_8dyn: Calibrated surrogate Modelo 8 Dynamic
  Atmospheric Coupling:
    - Renaudo (2022, 2023) evaporation parameter k, d_min(H), t_dep, d_dep
    - Chiu (Sirignano 2010) Group Evaporation Number G

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
println("CLICAP IN-FLIGHT DROPLET EVAPORATION: ALL MODELS & CASES")
println("================================================================================")

# 1. MATRIX OF 9 ENVIRONMENTAL CASES
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

# Export Table 1
df_env = DataFrame(
    Caso = [c.id for c in CLICAP_CASES],
    T_K = [c.T_K for c in CLICAP_CASES],
    T_C = [c.T_C for c in CLICAP_CASES],
    HR_pct = [c.HR for c in CLICAP_CASES],
    T_wb_C = [round(stull_twb(c.T_C, c.HR), digits=2) for c in CLICAP_CASES],
    psat_Pa = [round(psat_water(c.T_K), digits=1) for c in CLICAP_CASES],
    k_Renaudo_sm2 = [round(renaudo_evaporation_parameter(c.T_C, c.HR), digits=0) for c in CLICAP_CASES],
    d_min_50cm_um = [round(renaudo_extinction_diameter(0.5, c.T_C, c.HR)*1e6, digits=1) for c in CLICAP_CASES],
    d_min_100cm_um = [round(renaudo_extinction_diameter(1.0, c.T_C, c.HR)*1e6, digits=1) for c in CLICAP_CASES]
)

t1_csv = joinpath(DATA_DIR, "table1_environmental_cases.csv")
CSV.write(t1_csv, df_env)
println("Saved Table 1: $t1_csv")

# 2. SIMULATE ALL IN-FLIGHT COMBINATIONS
const DIAMETERS_UM = [50.0, 200.0, 600.0, 1200.0]
const P_ATM = 101325.0

const MODEL_LIST = [
    (key = "a", id = :d2_classic, profile = :uniforme, desc = "a (D² Clásica, Sh=2)"),
    (key = "b", id = :d2_stefan, profile = :uniforme, desc = "b (D² con Stefan)"),
    (key = "c", id = :d2_convective, profile = :uniforme, desc = "c (D² con Convección)"),
    (key = "d", id = :abramzon_sirignano, profile = :parabolico, desc = "d (A-S Parabólico)"),
    (key = "d_uni", id = :abramzon_sirignano, profile = :uniforme, desc = "d_uni (A-S Uniforme)"),
    (key = "d_zhou", id = :abramzon_sirignano, profile = :zhou_cubic, desc = "d_zhou (A-S Cúbico Zhou)"),
    (key = "revised", id = :d2_revised_dallabarba, profile = :uniforme, desc = "revised (D² Dalla Barba c/ Twb ventilado)"),
    (key = "vlasov", id = :vlasov_kinetic, profile = :uniforme, desc = "vlasov (Vlasov 2021 Cinético-Difusivo)"),
    (key = "surr_6b", id = :surrogate_6b, profile = :uniforme, desc = "surr_6b (Surrogante 6b Refit)"),
    (key = "surr_7", id = :surrogate_7, profile = :uniforme, desc = "surr_7 (Surrogante 7 Multi-Var)")
]

all_csv = joinpath(DATA_DIR, "comparison_inflight_all_cases.csv")
records_all = []
println("\nSimulating in-flight droplets across 9 cases × 4 diameters × $(length(MODEL_LIST)) models...")

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

        k_Renaudo = renaudo_evaporation_parameter(T_C, HR_pct)
        d_min_50 = renaudo_extinction_diameter(0.5, T_C, HR_pct)

        for d_um in DIAMETERS_UM
            d_m = d_um * 1e-6
            t_dep_50 = renaudo_deposition_time(d_m, 0.5, T_C, HR_pct)
            d_land_50 = renaudo_landing_diameter(d_m, 0.5, T_C, HR_pct)

            G_chiu = chiu_group_evaporation_number(100_000, d_m / 2.0, 0.05)
            regime_chiu = chiu_evaporation_regime(G_chiu)

            for m_info in MODEL_LIST
                m_key = m_info.key
                m_id = m_info.id
                m_prof = m_info.profile
                m_desc = m_info.desc

                t_evap = 0.0
                K_t0_mm2_s = 0.0

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

                    alg0 = calc_droplet_algebraic(m0, T_init_C, T_init_C, params)
                    K_t0_mm2_s = (4.0 * alg0.m_dot / (π * alg0.rho_l * alg0.d)) * 1e6

                    t_est = (d_um / 50.0)^2 * 10.0
                    tspan = (0.0, max(50.0, t_est * 4.0))

                    sol = solve_droplet(params, [m0, T_init_C], tspan)
                    t_evap = sol.t[end]

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
                    K_t0_mm2_s = (d_m * 1e3)^2 / t_evap

                elseif m_id in (:surrogate_6b, :surrogate_7)
                    stype = m_id == :surrogate_6b ? :modelo_6b_refit : :modelo_7
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
                    K_t0_mm2_s = (d_um * 1e-3)^2 / t_evap
                end

                K_avg_mm2_s = (d_um * 1e-3)^2 / t_evap

                push!(records_all, (
                    Caso = case.id,
                    T_C = T_C,
                    HR_pct = HR_pct,
                    d0_um = d_um,
                    Modelo_Key = m_key,
                    Modelo_Desc = m_desc,
                    t_evap_s = round(t_evap, digits=3),
                    K_t0_mm2_s = round(K_t0_mm2_s, digits=6),
                    K_avg_mm2_s = round(K_avg_mm2_s, digits=6),
                    k_Renaudo = round(k_Renaudo, digits=0),
                    d_min_50cm_um = round(d_min_50 * 1e6, digits=2),
                    t_dep_50cm_s = isinf(t_dep_50) ? NaN : round(t_dep_50, digits=3),
                    d_land_50cm_um = round(d_land_50 * 1e6, digits=2),
                    G_chiu = round(G_chiu, digits=3),
                    Regimen_Chiu = string(regime_chiu)
                ))
            end
        end
    end

    df_all = DataFrame(records_all)
    CSV.write(all_csv, df_all)
    println("Saved full dataset: $all_csv ($(nrow(df_all)) rows)")

# 3. DIRECT PAPER TABLE 2 REPRODUCTION
# Filters Case 1 and Case 7 for models a, b, c, d
sub_t2 = filter(r -> r.Caso in (1, 7) && r.Modelo_Key in ("a", "b", "c", "d"), df_all)

# Pivot to match paper table format
df_pivoted_t2 = unstack(sub_t2, [:Modelo_Key, :Modelo_Desc, :d0_um], :Caso, :t_evap_s, renamecols=c -> Symbol("t_Caso$(c)_Julia_s"))

# Reference values from paper Table 2
paper_ref = Dict(
    ("a", 50.0) => (1.26, 0.18),
    ("b", 50.0) => (1.27, 0.29),
    ("c", 50.0) => (1.19, 0.17),
    ("d", 50.0) => (3.13, 1.35),
    ("a", 200.0) => (20.4, 2.95),
    ("b", 200.0) => (20.4, 2.91),
    ("c", 200.0) => (14.4, 2.15),
    ("d", 200.0) => (36.1, 15.6),
    ("a", 600.0) => (184.0, 27.0),
    ("b", 600.0) => (183.0, 26.1),
    ("c", 600.0) => (78.0, 11.8),
    ("d", 600.0) => (190.0, 83.0),
    ("a", 1200.0) => (735.0, 106.0),
    ("b", 1200.0) => (733.0, 105.0),
    ("c", 1200.0) => (203.0, 31.0),
    ("d", 1200.0) => (488.0, 214.0)
)

df_pivoted_t2.t_Caso1_Paper_s = [paper_ref[(r.Modelo_Key, r.d0_um)][1] for r in eachrow(df_pivoted_t2)]
df_pivoted_t2.t_Caso7_Paper_s = [paper_ref[(r.Modelo_Key, r.d0_um)][2] for r in eachrow(df_pivoted_t2)]

t2_csv = joinpath(DATA_DIR, "table2_inflight_cases1_and_7.csv")
CSV.write(t2_csv, df_pivoted_t2)
println("Saved Table 2 Comparison: $t2_csv")

# 4. PLOTTING FIGURES
println("\nGenerating publication figures...")

# Figure A: In-flight lifetimes vs diameter for Case 1 (10°C, 20%) and Case 7 (40°C, 20%)
p_life = plot(
    title = "In-Flight Droplet Lifetimes: Model Comparison (Cases 1 & 7)",
    xlabel = "Initial Droplet Diameter d₀ [µm]",
    ylabel = "Evaporation Lifetime t_evap [s]",
    xscale = :log10,
    yscale = :log10,
    legend = :topleft,
    framestyle = :box,
    grid = true,
    dpi = 300,
    size = (800, 600)
)

colors_dict = Dict(
    "a" => :royalblue,
    "b" => :cyan4,
    "c" => :forestgreen,
    "d" => :crimson,
    "d_zhou" => :magenta,
    "revised" => :darkorange,
    "surr_6b" => :purple
)

for m_k in ["a", "c", "d", "d_zhou", "revised"]
    sub_c1 = filter(r -> r.Caso == 1 && r.Modelo_Key == m_k, df_all)
    sub_c7 = filter(r -> r.Caso == 7 && r.Modelo_Key == m_k, df_all)

    lbl = sub_c1.Modelo_Desc[1]
    clr = colors_dict[m_k]

    plot!(p_life, sub_c1.d0_um, sub_c1.t_evap_s, label="$lbl (C1 10°C)", color=clr, lw=2, marker=:circle)
    plot!(p_life, sub_c7.d0_um, sub_c7.t_evap_s, label="$lbl (C7 40°C)", color=clr, lw=2, ls=:dash, marker=:square)
end

fig_a_png = joinpath(FIG_DIR, "fig_inflight_lifetimes_by_diameter.png")
savefig(p_life, fig_a_png)
println("Saved Figure A: $fig_a_png")

# Figure B: Initial Evaporation Rate K(t=0) vs Relative Humidity (Cases 4, 5, 6: T=25°C, d0=200 µm)
sub_f3 = filter(r -> r.Caso in (4, 5, 6) && r.d0_um == 200.0 && r.Modelo_Key in ("a", "b", "c", "d", "d_zhou", "revised"), df_all)

p_k0 = plot(
    title = "Initial Evaporation Rate K(t=0) vs RH (d₀ = 200 µm, T = 25 °C)",
    xlabel = "Relative Humidity HR [%]",
    ylabel = "K(t=0) [mm²/s]",
    legend = :topright,
    framestyle = :box,
    grid = true,
    dpi = 300,
    size = (750, 500)
)

for m_k in ["a", "b", "c", "d", "d_zhou", "revised"]
    sub_m = filter(r -> r.Modelo_Key == m_k, sub_f3)
    lbl = sub_m.Modelo_Desc[1]
    clr = get(colors_dict, m_k, :black)
    plot!(p_k0, sub_m.HR_pct, sub_m.K_t0_mm2_s, label=lbl, color=clr, lw=2, marker=:circle)
end

fig_b_png = joinpath(FIG_DIR, "fig_inflight_K_initial_vs_HR.png")
savefig(p_k0, fig_b_png)
println("Saved Figure B: $fig_b_png")

# Figure C: Renaudo In-Flight Extinction Cut-Off Diameter d_min vs Weather (Cases 1 to 9)
p_ext = plot(
    title = "Cut-Off Extinction Diameter d_min vs RH across Temperatures",
    xlabel = "Relative Humidity HR [%]",
    ylabel = "Cut-Off Diameter d_min [µm] (H = 0.5 m)",
    legend = :topright,
    framestyle = :box,
    grid = true,
    dpi = 300,
    size = (750, 500)
)

for t_c in [10.0, 25.0, 40.0]
    sub_t = filter(r -> r.T_C == t_c, df_env)
    lbl = @sprintf("T = %.0f °C", t_c)
    plot!(p_ext, sub_t.HR_pct, sub_t.d_min_50cm_um, label=lbl, lw=2.5, marker=:diamond, ms=6)
end

fig_c_png = joinpath(FIG_DIR, "fig_inflight_extinction_dmin_vs_weather.png")
savefig(p_ext, fig_c_png)
println("Saved Figure C: $fig_c_png")

println("\nIn-flight script execution completed successfully!")
