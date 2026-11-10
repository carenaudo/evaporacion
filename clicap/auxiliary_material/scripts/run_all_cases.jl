"""
Master Reproduction & Execution Script - Auxiliary Material
Executes both in-flight and sessile droplet case studies across all models and scenarios,
synthesizes comparative cross-domain metrics, and generates publication artifacts.

Usage:
  julia --project=julia clicap/auxiliary_material/scripts/run_all_cases.jl
"""

ENV["GKSwstype"] = "100"
println("================================================================================")
println("STARTING COMPLETE CLICAP REPRODUCTION & MODEL COMPARISON SUITE")
println("================================================================================")

const SCRIPT_DIR = @__DIR__

# Load Evaporacion.jl package once at master level
if !isdefined(Main, :Evaporacion)
    include(joinpath(SCRIPT_DIR, "..", "..", "..", "julia", "src", "Evaporacion.jl"))
    using .Evaporacion
end

# 1. Execute In-Flight Droplet Simulations
println("\n>>> [1/2] RUNNING IN-FLIGHT SIMULATIONS...")
include(joinpath(SCRIPT_DIR, "run_inflight_cases.jl"))

# 2. Execute Sessile Droplet Simulations
println("\n>>> [2/2] RUNNING SESSILE SIMULATIONS...")
include(joinpath(SCRIPT_DIR, "run_sessile_cases.jl"))

# 3. Cross-Domain Comparative Synthesis: Flying vs Sessile for d0 = 200 µm
println("\n>>> GENERATING CROSS-DOMAIN SYNTHESIS (FLYING VS SESSILE)...")
using DataFrames
using CSV
using Plots
using Printf

const DATA_DIR = normpath(joinpath(SCRIPT_DIR, "..", "data"))
const FIG_DIR  = normpath(joinpath(SCRIPT_DIR, "..", "figures"))

df_inf = CSV.read(joinpath(DATA_DIR, "comparison_inflight_all_cases.csv"), DataFrame)
df_ses = CSV.read(joinpath(DATA_DIR, "comparison_sessile_all_cases.csv"), DataFrame)

# For 200 µm droplets across all 9 cases
records_cross = []

for c_id in 1:9
    # In-flight representative models
    sub_inf_c = filter(r -> r.Caso == c_id && r.d0_um == 200.0, df_inf)
    t_d2_classic = filter(r -> r.Modelo_Key == "a", sub_inf_c).t_evap_s[1]
    t_d2_conv    = filter(r -> r.Modelo_Key == "c", sub_inf_c).t_evap_s[1]
    t_as_para    = filter(r -> r.Modelo_Key == "d", sub_inf_c).t_evap_s[1]
    t_d_zhou     = filter(r -> r.Modelo_Key == "d_zhou", sub_inf_c).t_evap_s[1]
    t_vlasov     = filter(r -> r.Modelo_Key == "vlasov", sub_inf_c).t_evap_s[1]
    d_min_50     = sub_inf_c.d_min_50cm_um[1]
    t_dep_50     = sub_inf_c.t_dep_50cm_s[1]
    d_land_50    = sub_inf_c.d_land_50cm_um[1]

    # Sessile representative models at theta0 = 40°
    sub_ses_c = filter(r -> r.Caso == c_id && r.theta0_deg == 40.0, df_ses)
    t_ca = sub_ses_c.t_AC_s[1]
    t_cr = sub_ses_c.t_RC_s[1]
    t_ss = sub_ses_c.t_SS_s[1]
    t_sj = sub_ses_c.t_SJ_s[1]
    t_cool = sub_ses_c.t_Cooled_s[1]

    push!(records_cross, (
        Caso = c_id,
        T_C = sub_inf_c.T_C[1],
        HR_pct = sub_inf_c.HR_pct[1],
        d0_um = 200.0,
        t_Vuelo_D2_Clasica_s = t_d2_classic,
        t_Vuelo_Convectiva_s = t_d2_conv,
        t_Vuelo_Abramzon_Sirignano_s = t_as_para,
        t_Vuelo_Zhou_Cubico_s = t_d_zhou,
        t_Vuelo_Vlasov_Knudsen_s = t_vlasov,
        d_min_Extincion_50cm_um = d_min_50,
        t_dep_Sedimentacion_50cm_s = t_dep_50,
        d_land_Suelo_50cm_um = d_land_50,
        t_Sesil_AC_40deg_s = t_ca,
        t_Sesil_RC_40deg_s = t_cr,
        t_Sesil_SS_40deg_s = t_ss,
        t_Sesil_SJ_40deg_s = t_sj,
        t_Sesil_Enfriada_40deg_s = t_cool,
        Ratio_VueloAS_vs_SesilRC = round(t_as_para / t_cr, digits=2)
    ))
end

df_cross = DataFrame(records_cross)
cross_csv = joinpath(DATA_DIR, "comparison_cross_domain_flying_vs_sessile.csv")
CSV.write(cross_csv, df_cross)
println("Saved Cross-Domain Dataset: $cross_csv")

# Generate Cross-Domain Figure
p_cross = plot(
    title = "Flying vs Sessile Lifetimes (d₀ = 200 µm, θ₀ = 40°)",
    xlabel = "CLICAP Environmental Case Index (1 to 9)",
    ylabel = "Evaporation Lifetime [s]",
    yscale = :log10,
    legend = :topright,
    framestyle = :box,
    grid = true,
    dpi = 300,
    size = (850, 520)
)

cases_x = 1:9
plot!(p_cross, cases_x, df_cross.t_Vuelo_D2_Clasica_s, label="In-Flight D² Clásica (Estanco)", lw=2, color=:royalblue, marker=:circle)
plot!(p_cross, cases_x, df_cross.t_Vuelo_Abramzon_Sirignano_s, label="In-Flight Abramzon-Sirignano", lw=2.5, color=:crimson, marker=:square)
plot!(p_cross, cases_x, df_cross.t_Sesil_RC_40deg_s, label="Sessile RC (Pinned θ₀=40°)", lw=2.5, color=:forestgreen, marker=:diamond)
plot!(p_cross, cases_x, df_cross.t_Sesil_AC_40deg_s, label="Sessile AC (Receding θ₀=40°)", lw=2, color=:darkorange, marker=:triangle)
plot!(p_cross, cases_x, df_cross.t_Sesil_Enfriada_40deg_s, label="Sessile RC Substrate Cooled", lw=2, ls=:dash, color=:purple)

case_labels = [@sprintf("C%d\n%.0f°C\n%.0f%%", r.Caso, r.T_C, r.HR_pct) for r in eachrow(df_cross)]
xticks!(p_cross, cases_x, case_labels)

fig_cross_png = joinpath(FIG_DIR, "fig_flying_vs_sessile_all_cases.png")
savefig(p_cross, fig_cross_png)
println("Saved Cross-Domain Figure: $fig_cross_png")

println("\n================================================================================")
println("ALL SIMULATION CASES & COMPARATIVE DATASETS SUCCESSFULLY GENERATED!")
println("Output Data Directory   : $DATA_DIR")
println("Output Figures Directory: $FIG_DIR")
ENV["GKSwstype"] = "100"
println("================================================================================")
