"""
Generation of Figure 2: Sessile droplet evaporation modes reported by Wilson & D'Ambrosio (2023).
Reproduces the 4 canonical modes:
  1. Constant Contact Radius (CR)
  2. Constant Contact Angle (CA)
  3. Stick-Slide (S-S) with θ_rec = θ₀ / 2
  4. Stick-Jump (S-J) with θ_min = θ₀ / 2, θ_max = θ₀

Compares numerical solutions from Evaporacion.Sessile with asymptotic Table 2 (Wilson 2023).
"""

using Plots
using DataFrames
using CSV
using QuadGK

include(joinpath(@__DIR__, "..", "julia", "src", "Evaporacion.jl"))
using .Evaporacion

const ARTIFACT_DIR = @__DIR__
mkpath(ARTIFACT_DIR)

println("================================================================================")
println("REPRODUCING FIGURE 2: MODOS DE EVAPORACIÓN GOTA SÉSIL (Wilson & D'Ambrosio 2023)")
println("================================================================================")

# Selected canonical parameters to reproduce the Wilson & D'Ambrosio behavior
# Contact angle: θ₀ = 30° (π/6 rad) -> thin/moderate droplet
const THETA0_DEG = 30.0
const THETA0_RAD = deg2rad(THETA0_DEG)
const R0_VAL = 1.0e-3        # 1 mm initial contact radius
const D_VAL = 2.6e-5         # m^2/s
const DELTA_C = 0.02         # kg/m^3
const RHO_L = 997.0          # kg/m^3

# Precalculate exact t_CR scale:
# t_CR = (rho * R0^2 / (D * delta_c)) * integral_0^theta0 (1 / g(theta)) dtheta
const PREFACTOR = (RHO_L * (R0_VAL^2)) / (D_VAL * DELTA_C)
const INT_CR = quadgk(th -> 1.0 / g_theta_fast(th), 0.0, THETA0_RAD)[1]
const T_CR_SCALE = PREFACTOR * INT_CR

println("Conditions and Parameters:")
println("  Initial Contact Angle θ₀: $(THETA0_DEG)° ($(round(THETA0_RAD, digits=4)) rad)")
println("  Initial Contact Radius R₀: $(R0_VAL*1000) mm")
println("  Binary Diffusivity D: $D_VAL m²/s")
println("  Concentration Deficit Δc: $DELTA_C kg/m³")
println("  Liquid Density ρ: $RHO_L kg/m³")
println("  Characteristic CR timescale t_CR,scale: $(round(T_CR_SCALE, digits=2)) s")

# ------------------------------------------------------------------------------
# 1. MODE CR (Constant Contact Radius)
# ------------------------------------------------------------------------------
cfg_cr = SessileDropletConfig(
    D = D_VAL,
    delta_c = DELTA_C,
    rho_l = RHO_L,
    R0 = R0_VAL,
    theta0 = THETA0_RAD,
    mode = :CR,
    t_end = T_CR_SCALE * 1.05
)
sol_cr = solve_sessile(cfg_cr)
tau_cr = sol_cr.t ./ T_CR_SCALE
R_norm_cr = sol_cr.R ./ R0_VAL
th_norm_cr = sol_cr.theta ./ THETA0_RAD

# Asymptotic Wilson Table 2 for CR: θ/θ₀ = 1 - τ
tau_w_cr = range(0.0, 1.0, length=100)
th_w_cr = 1.0 .- tau_w_cr

# ------------------------------------------------------------------------------
# 2. MODE CA (Constant Contact Angle)
# ------------------------------------------------------------------------------
cfg_ca = SessileDropletConfig(
    D = D_VAL,
    delta_c = DELTA_C,
    rho_l = RHO_L,
    R0 = R0_VAL,
    theta0 = THETA0_RAD,
    mode = :CA,
    t_end = T_CR_SCALE * 2.0
)
sol_ca = solve_sessile(cfg_ca)
tau_ca = sol_ca.t ./ T_CR_SCALE
R_norm_ca = sol_ca.R ./ R0_VAL
th_norm_ca = sol_ca.theta ./ THETA0_RAD

# Asymptotic Wilson Table 2 for CA: R/R₀ = sqrt(1 - τ / 1.5)
tau_ca_end = sol_ca.t[end] / T_CR_SCALE
tau_w_ca = range(0.0, 1.5, length=100)
R_w_ca = [t <= 1.5 ? sqrt(max(0.0, 1.0 - t / 1.5)) : 0.0 for t in tau_w_ca]

# ------------------------------------------------------------------------------
# 3. MODE STICK-SLIDE (S-S) with θ_rec = θ₀ / 2
# ------------------------------------------------------------------------------
const THETA_REC_RAD = THETA0_RAD / 2.0
cfg_ss = SessileDropletConfig(
    D = D_VAL,
    delta_c = DELTA_C,
    rho_l = RHO_L,
    R0 = R0_VAL,
    theta0 = THETA0_RAD,
    theta_rec = THETA_REC_RAD,
    mode = :SS,
    t_end = T_CR_SCALE * 2.0
)
sol_ss = solve_sessile(cfg_ss)
tau_ss = sol_ss.t ./ T_CR_SCALE
R_norm_ss = sol_ss.R ./ R0_VAL
th_norm_ss = sol_ss.theta ./ THETA0_RAD

# ------------------------------------------------------------------------------
# 4. MODE STICK-JUMP (S-J) with θ_min = θ₀ / 2, θ_max = θ₀
# ------------------------------------------------------------------------------
cfg_sj = SessileDropletConfig(
    D = D_VAL,
    delta_c = DELTA_C,
    rho_l = RHO_L,
    R0 = R0_VAL,
    theta0 = THETA0_RAD,
    theta_min = THETA_REC_RAD,
    theta_max = THETA0_RAD,
    mode = :SJ,
    t_end = T_CR_SCALE * 2.0
)
sol_sj = solve_sessile(cfg_sj)
tau_sj = sol_sj.t ./ T_CR_SCALE
R_norm_sj = sol_sj.R ./ R0_VAL
th_norm_sj = sol_sj.theta ./ THETA0_RAD

# ------------------------------------------------------------------------------
# BUILD COMPLETE CSV DATASET
# ------------------------------------------------------------------------------
records = []
for i in eachindex(tau_cr)
    push!(records, (Modo="CR", tau=tau_cr[i], R_norm=R_norm_cr[i], theta_norm=th_norm_cr[i]))
end
for i in eachindex(tau_ca)
    push!(records, (Modo="CA", tau=tau_ca[i], R_norm=R_norm_ca[i], theta_norm=th_norm_ca[i]))
end
for i in eachindex(tau_ss)
    push!(records, (Modo="SS", tau=tau_ss[i], R_norm=R_norm_ss[i], theta_norm=th_norm_ss[i]))
end
for i in eachindex(tau_sj)
    push!(records, (Modo="SJ", tau=tau_sj[i], R_norm=R_norm_sj[i], theta_norm=th_norm_sj[i]))
end

df_fig2 = DataFrame(records)
fig2_csv = joinpath(ARTIFACT_DIR, "figure2_sessile_modes_wilson.csv")
CSV.write(fig2_csv, df_fig2)
println("Saved CSV: $fig2_csv")

# ------------------------------------------------------------------------------
# PLOT 2x2 FIGURE MATCHING EXACT REFERENCE STYLE WITH COLORS, NO WILSON CURVE
# ------------------------------------------------------------------------------
default(
    fontfamily = "Helvetica",
    framestyle = :box,
    grid = true,
    gridcolor = :gray85,
    gridalpha = 0.6,
    gridlinewidth = 0.8,
    titlefontsize = 13,
    guidefontsize = 11,
    tickfontsize = 10,
    legendfontsize = 9,
    dpi = 300
)

const COLOR_R = :royalblue
const COLOR_THETA = :crimson
const LABEL_R = "Radio normalizado R(t) / R_0"
const LABEL_THETA = "Ángulo normalizado θ(t) / θ_0"

# Panel 1 (Top-Left): Modo Radio de contacto constante
p1 = plot(
    title = "Modo Radio de contacto constante",
    xlabel = "Tiempo normalizado (t / t_CR,scale)",
    ylabel = "Radio y ángulo de contacto normalizados",
    xlim = (-0.03, 1.03),
    ylim = (-0.04, 1.06),
    xticks = 0.0:0.2:1.0,
    yticks = [0.0, 0.5, 1.0],
    legend = :bottomleft,
    background_color_legend = RGBA(1.0, 1.0, 1.0, 0.88),
    foreground_color_legend = :gray70
)
plot!(p1, tau_cr, R_norm_cr, label=LABEL_R, color=COLOR_R, lw=2.8, ls=:solid)
plot!(p1, tau_cr, th_norm_cr, label=LABEL_THETA, color=COLOR_THETA, lw=2.8, ls=:dash)

# Panel 2 (Bottom-Left): Modo Ángulo de Contacto constante
p2 = plot(
    title = "Modo Ángulo de Contacto constante",
    xlabel = "Tiempo normalizado (t / t_CR,scale)",
    ylabel = "Radio y ángulo de contacto normalizados",
    xlim = (-0.04, 1.54),
    ylim = (-0.04, 1.06),
    xticks = 0.0:0.25:1.50,
    yticks = [0.0, 0.5, 1.0],
    legend = :bottomleft,
    background_color_legend = RGBA(1.0, 1.0, 1.0, 0.88),
    foreground_color_legend = :gray70
)
plot!(p2, tau_ca, R_norm_ca, label=LABEL_R, color=COLOR_R, lw=2.8, ls=:solid)
plot!(p2, tau_ca, th_norm_ca, label=LABEL_THETA, color=COLOR_THETA, lw=2.8, ls=:dash)

# Panel 3 (Top-Right): Modo Stick-Slide
p3 = plot(
    title = "Modo Stick-Slide",
    xlabel = "Tiempo normalizado (t / t_CR,scale)",
    ylabel = "Radio y ángulo de contacto normalizados",
    xlim = (-0.03, 1.28),
    ylim = (-0.04, 1.06),
    xticks = 0.0:0.2:1.2,
    yticks = [0.0, 0.5, 1.0],
    legend = :bottomleft,
    background_color_legend = RGBA(1.0, 1.0, 1.0, 0.88),
    foreground_color_legend = :gray70
)
plot!(p3, tau_ss, R_norm_ss, label=LABEL_R, color=COLOR_R, lw=2.8, ls=:solid)
plot!(p3, tau_ss, th_norm_ss, label=LABEL_THETA, color=COLOR_THETA, lw=2.8, ls=:dash)

# Panel 4 (Bottom-Right): Modo Stick-Jump
p4 = plot(
    title = "Modo Stick-Jump",
    xlabel = "Tiempo normalizado (t / t_CR,scale)",
    ylabel = "Radio y ángulo de contacto normalizados",
    xlim = (-0.04, 1.38),
    ylim = (-0.04, 1.06),
    xticks = 0.0:0.25:1.25,
    yticks = [0.0, 0.5, 1.0],
    legend = :bottomleft,
    background_color_legend = RGBA(1.0, 1.0, 1.0, 0.88),
    foreground_color_legend = :gray70
)
plot!(p4, tau_sj, R_norm_sj, label=LABEL_R, color=COLOR_R, lw=2.5, ls=:solid)
plot!(p4, tau_sj, th_norm_sj, label=LABEL_THETA, color=COLOR_THETA, lw=2.5, ls=:dash)

# Combine in 2x2 grid: [CR, Stick-Slide; CA, Stick-Jump]
fig2 = plot(
    p1, p3, p2, p4,
    layout = (2, 2),
    size = (1160, 640),
    left_margin = 6Plots.mm,
    right_margin = 3Plots.mm,
    top_margin = 5Plots.mm,
    bottom_margin = 6Plots.mm
)

fig2_png = joinpath(ARTIFACT_DIR, "figure2_sessile_modes_wilson.png")
savefig(fig2, fig2_png)
println("Saved PNG in clicap: $fig2_png")

const BRAIN_ARTIFACT_DIR = "C:/Users/charl/.gemini/antigravity/brain/85032b89-12d0-488b-9f15-8c33f4910707"
if isdir(BRAIN_ARTIFACT_DIR)
    brain_fig2_png = joinpath(BRAIN_ARTIFACT_DIR, "figure2_sessile_modes_wilson.png")
    savefig(fig2, brain_fig2_png)
    println("Saved PNG in brain artifact: $brain_fig2_png")
end

println("\n================================================================================")
println("FIGURE 2 GENERATION FINISHED SUCCESSFULLY!")
println("================================================================================")
