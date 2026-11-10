"""
Case study simulation and figure generation script.
Runs all models in Julia for:
  T_inf = 35.0 °C
  HR    = 0.40 (40%)
  P_atm = 101325 Pa
  U_rel = 0.5 m/s
  d0    = 1.0 mm (V0 ≈ 0.5236 mm³)
  T0    = 22.0 °C
"""

using Plots
using DataFrames
using QuadGK

# Load Evaporacion suite
include(joinpath(@__DIR__, "..", "julia", "src", "Evaporacion.jl"))
using .Evaporacion

# Define artifact directory
const ARTIFACT_DIR = @__DIR__

# Case parameters
const T_INF = 35.0
const HR_FRAC = 0.40
const P_ATM_VAL = 101325.0
const U_REL = 0.5
const D0_MM = 1.0
const D0_M = D0_MM * 1e-3
const T0_C = 22.0

println("==================================================")
println("RUNNING SINGLE CASE STUDY (JULIA)")
println("Ambient: T_inf = $T_INF °C, HR = $(HR_FRAC*100)%, P = $P_ATM_VAL Pa, U = $U_REL m/s")
println("Droplet: d0 = $D0_MM mm, T0 = $T0_C °C")
println("==================================================")

# ----------------------------------------------------
# 1. RIGOROUS FLYING DROPLET
# ----------------------------------------------------
println("\n1. Solving Rigorous Flying Droplet...")
T0_K = T0_C + Evaporacion.Common.Constants.T_ZERO_C
rho_l0 = liquid_density_water(T0_K)
V0_m3 = (π / 6.0) * (D0_M^3)
m0 = rho_l0 * V0_m3

# Far-field vapor mass fraction
T_inf_K = T_INF + Evaporacion.Common.Constants.T_ZERO_C
psat_inf = psat_water(T_inf_K)
p_v_inf = HR_FRAC * psat_inf
p_a_inf = P_ATM_VAL - p_v_inf
Y_inf = (p_v_inf * Evaporacion.Common.Constants.M_WATER) / 
        (p_v_inf * Evaporacion.Common.Constants.M_WATER + p_a_inf * Evaporacion.Common.Constants.M_AIR)

params_rig = DropletParams(
    T_inf = T_INF,
    U_rel = U_REL,
    Y_inf = Y_inf,
    P_atm = P_ATM_VAL,
    profile = :parabolico,
    correlation = :frossling
)

sol_rig = solve_droplet(params_rig, [m0, T0_C], (0.0, 800.0))
df_rig = postprocess_trajectory(sol_rig, params_rig)
t_evap_rig = df_rig.t[end]
println("   Rigorous evaporation lifetime: $(round(t_evap_rig, digits=2)) s (final d: $(round(df_rig.d[end]*1000, digits=4)) mm)")

# ----------------------------------------------------
# 2. CALIBRATED SURROGATES
# ----------------------------------------------------
println("\n2. Solving Calibrated Surrogates...")
T_wb_val = stull_twb(T_INF, HR_FRAC * 100.0)
println("   Stull Wet-Bulb Temperature: $(round(T_wb_val, digits=2)) °C")

# Modelo 6b Refit
cfg_6b = SurrogateConfig(
    model_type = :modelo_6b_refit,
    d0_mm = D0_MM,
    T_inf = T_INF,
    T_wb = T_wb_val,
    U = U_REL,
    HR = HR_FRAC,
    T0 = T0_C,
    t_end = 800.0
)
sol_6b = solve_surrogate(cfg_6b)
t_evap_6b = sol_6b.t[end]
println("   Modelo 6b Refit lifetime: $(round(t_evap_6b, digits=2)) s (final V: $(round(sol_6b.u[end][1], digits=6)) mm³)")

# Modelo 8 Dynamic (coupled V & T)
cfg_m8 = SurrogateConfig(
    model_type = :modelo_8_dynamic,
    d0_mm = D0_MM,
    T_inf = T_INF,
    T_wb = T_wb_val,
    U = U_REL,
    HR = HR_FRAC,
    T0 = T0_C,
    t_end = 800.0
)
sol_m8 = solve_surrogate(cfg_m8)
t_evap_m8 = sol_m8.t[end]
println("   Modelo 8 Dynamic lifetime: $(round(t_evap_m8, digits=2)) s (final V: $(round(sol_m8.u[end][1], digits=6)) mm³)")

# ----------------------------------------------------
# 3. 1D RADIAL CONDUCTION MOVING-BOUNDARY PDE
# ----------------------------------------------------
println("\n3. Solving 1D Radial Moving-Boundary PDE...")
cfg_pde = RadialPDEConfig(
    n_nodes = 31,
    d0 = D0_M,
    T_drop_0 = T0_C,
    T_air = T_INF,
    RH = HR_FRAC,
    with_surfactants = false,
    with_marangoni = true,
    with_natural_convection = true,
    t_end = 800.0
)
sol_pde = solve_radial_pde(cfg_pde)
t_evap_pde = sol_pde.t[end]
println("   1D Radial PDE lifetime: $(round(t_evap_pde, digits=2)) s (final d: $(round(sol_pde.R[end]*2000, digits=4)) mm)")

# ----------------------------------------------------
# 4. SESSILE DROPLET (Wilson & D'Ambrosio 2023 / Popov 2005)
# ----------------------------------------------------
println("\n4. Solving Sessile Droplets (CR, CA, SS, Cooled)...")
# Ambient vapor concentrations
c_sat = (psat_inf * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * T_inf_K)
c_inf = HR_FRAC * c_sat
delta_c_val = c_sat - c_inf
D_val = gas_diffusivity_water_air(T_inf_K, P_ATM_VAL)

# For θ0 = 90° (π/2), droplet is a hemisphere:
# V0 = (2/3)*π*R0^3. For same volume V0_m3: R0 = (3*V0/(2*π))^(1/3)
theta0_rad = 0.5 * π
R0_sessile = (3.0 * V0_m3 / (2.0 * π))^(1.0 / 3.0)
println("   Sessile equivalent hemisphere contact radius R0: $(round(R0_sessile*1000, digits=3)) mm")

# Constant Radius (CR)
cfg_cr = SessileDropletConfig(
    D = D_val,
    delta_c = delta_c_val,
    rho_l = rho_l0,
    R0 = R0_sessile,
    theta0 = theta0_rad,
    mode = :CR,
    t_end = 1500.0
)
sol_cr = solve_sessile(cfg_cr)
t_evap_cr = sol_cr.t[end]
println("   Sessile CR lifetime: $(round(t_evap_cr, digits=2)) s (exact analytical: $(round(sol_cr.t_evap_analytical, digits=2)) s)")

# Constant Angle (CA)
cfg_ca = SessileDropletConfig(
    D = D_val,
    delta_c = delta_c_val,
    rho_l = rho_l0,
    R0 = R0_sessile,
    theta0 = theta0_rad,
    mode = :CA,
    t_end = 1500.0
)
sol_ca = solve_sessile(cfg_ca)
t_evap_ca = sol_ca.t[end]
println("   Sessile CA lifetime: $(round(t_evap_ca, digits=2)) s (exact analytical: $(round(sol_ca.t_evap_analytical, digits=2)) s)")

# Stick-Slip (SS) with receding angle theta_rec = 30° (π/6)
cfg_ss = SessileDropletConfig(
    D = D_val,
    delta_c = delta_c_val,
    rho_l = rho_l0,
    R0 = R0_sessile,
    theta0 = theta0_rad,
    theta_rec = π / 6.0,
    mode = :SS,
    t_end = 1500.0
)
sol_ss = solve_sessile(cfg_ss)
t_evap_ss = sol_ss.t[end]
println("   Sessile SS lifetime: $(round(t_evap_ss, digits=2)) s (exact analytical: $(round(sol_ss.t_evap_analytical, digits=2)) s)")

# Substrate Cooling (CR with glass substrate: hs = 1 mm, ks = 1.4 W/(m K))
L_heat = latent_heat_water(T0_K)
dcsat_dT = (psat_water(T_inf_K + 0.1) - psat_water(T_inf_K - 0.1)) / 0.2 * (Evaporacion.Common.Constants.M_WATER / (Evaporacion.Common.Constants.R_GAS * T_inf_K))
thermal_glass = SessileThermalParams(
    L = L_heat,
    h_s = 1.0e-3,
    k_s_eff = 1.4,
    dcsat_dT = dcsat_dT
)
cfg_cooled = SessileDropletConfig(
    D = D_val,
    delta_c = delta_c_val,
    rho_l = rho_l0,
    R0 = R0_sessile,
    theta0 = theta0_rad,
    thermal = thermal_glass,
    mode = :CR,
    t_end = 1500.0
)
sol_cooled = solve_sessile(cfg_cooled)
t_evap_cooled = sol_cooled.t[end]
println("   Sessile Substrate Cooled lifetime: $(round(t_evap_cooled, digits=2)) s")

# ----------------------------------------------------
# 5. GENERATE PUBLICATION-QUALITY PLOTS
# ----------------------------------------------------
println("\n5. Generating Plots via Plots.jl...")
default(fontfamily="Helvetica", dpi=300, grid=true, framestyle=:box, lw=2)

# FIG 1: Flying Droplet d²(t) and Volume V(t)
V0_mm3 = (π / 6.0) * (D0_MM^3)
d2_0 = D0_MM^2

# Rigorous
t_rig = df_rig.t
d2_rig = (df_rig.d .* 1000.0).^2
V_rig = (π / 6.0) .* ((df_rig.d .* 1000.0).^3)

# 6b refit
t_6b = sol_6b.t
V_6b = [u[1] for u in sol_6b.u]
d2_6b = (6.0 .* V_6b ./ π).^(2.0 / 3.0)

# Modelo 8
t_m8 = sol_m8.t
V_m8 = [u[1] for u in sol_m8.u]
d2_m8 = (6.0 .* V_m8 ./ π).^(2.0 / 3.0)

# 1D Radial PDE
t_pde = sol_pde.t
R_pde_mm = sol_pde.R .* 1000.0
d2_pde = (2.0 .* R_pde_mm).^2
V_pde = (4.0 / 3.0) .* π .* (R_pde_mm.^3)

t_max_fly = max(t_evap_rig, t_evap_6b, t_evap_m8, t_evap_pde) * 1.05

p1a = plot(title="Flying Droplet d²-law Kinetics (T∞=35°C, HR=40%, U=0.5m/s)",
           xlabel="Time [s]", ylabel="d² / d₀² [-]", xlim=(0, t_max_fly), ylim=(0, 1.05))
plot!(p1a, t_rig, d2_rig ./ d2_0, label="Rigorous (Abramzon-Sirignano)", color=:navy, lw=2.5)
plot!(p1a, t_6b, d2_6b ./ d2_0, label="Modelo 6b Refit (Calibrated)", color=:crimson, lw=2, ls=:dash)
plot!(p1a, t_m8, d2_m8 ./ d2_0, label="Modelo 8 Dynamic (Coupled)", color=:darkgreen, lw=2, ls=:dot)
plot!(p1a, t_pde, d2_pde ./ d2_0, label="1D Radial Conduction PDE", color=:darkorange, lw=2, ls=:dashdot)

p1b = plot(xlabel="Time [s]", ylabel="Volume V(t) [mm³]", xlim=(0, t_max_fly), ylim=(0, V0_mm3*1.05))
plot!(p1b, t_rig, V_rig, label="Rigorous", color=:navy, lw=2.5)
plot!(p1b, t_6b, V_6b, label="Modelo 6b Refit", color=:crimson, lw=2, ls=:dash)
plot!(p1b, t_m8, V_m8, label="Modelo 8 Dynamic", color=:darkgreen, lw=2, ls=:dot)
plot!(p1b, t_pde, V_pde, label="1D Radial PDE", color=:darkorange, lw=2, ls=:dashdot)

fig1 = plot(p1a, p1b, layout=(2, 1), size=(850, 750))
fig1_path = joinpath(ARTIFACT_DIR, "fig1_flying_droplet_models.png")
savefig(fig1, fig1_path)
println("   Saved: $fig1_path")

# FIG 2: Droplet Temperature Evolution
p2 = plot(title="Flying Droplet Thermal Dynamics (T₀=22°C, T∞=35°C, HR=40%)",
          xlabel="Time [s]", ylabel="Temperature [°C]", xlim=(0, t_max_fly), ylim=(15, 36), legend=:bottomright)
plot!(p2, t_rig, df_rig.T_bar, label="Rigorous Mean T_bar(t)", color=:blue, lw=2.5)
plot!(p2, t_rig, df_rig.Ts, label="Rigorous Surface Ts(t)", color=:deepskyblue, lw=2, ls=:dash)
plot!(p2, t_m8, [u[2] for u in sol_m8.u], label="Modelo 8 Dynamic T(t)", color=:darkgreen, lw=2, ls=:dot)
# PDE center and surface
T_pde_center = sol_pde.T_profile[:, 1]
T_pde_surf = sol_pde.T_profile[:, end]
plot!(p2, t_pde, T_pde_center, label="1D PDE Center T(r=0, t)", color=:darkorange, lw=2)
plot!(p2, t_pde, T_pde_surf, label="1D PDE Surface T(r=R, t)", color=:gold, lw=2, ls=:dash)
hline!(p2, [T_wb_val], label="Stull Wet-Bulb Limit ($round(T_wb_val, digits=1)°C)", color=:gray, lw=1.5, ls=:dashdot)
hline!(p2, [T_INF], label="Ambient T∞ (35°C)", color=:red, lw=1.2, ls=:dot)

fig2_path = joinpath(ARTIFACT_DIR, "fig2_flying_droplet_temperatures.png")
savefig(p2, fig2_path)
println("   Saved: $fig2_path")

# FIG 3: Sessile Droplet Dynamics (CR, CA, SS Modes)
p3a = plot(title="Sessile Droplet Dynamics (Wilson & D'Ambrosio 2023 / Popov 2005)",
           xlabel="Time [s]", ylabel="Contact Radius R [mm]", xlim=(0, max(t_evap_cr, t_evap_ca)*1.05))
plot!(p3a, sol_cr.t, sol_cr.R .* 1000.0, label="Constant Radius (CR)", color=:navy, lw=2.5)
plot!(p3a, sol_ca.t, sol_ca.R .* 1000.0, label="Constant Angle (CA)", color=:darkgreen, lw=2.5, ls=:dash)
plot!(p3a, sol_ss.t, sol_ss.R .* 1000.0, label="Stick-Slip (SS, θ_rec=30°)", color=:purple, lw=2, ls=:dashdot)
plot!(p3a, sol_cooled.t, sol_cooled.R .* 1000.0, label="Substrate Cooled (CR, Glass)", color=:crimson, lw=2, ls=:dot)

p3b = plot(xlabel="Time [s]", ylabel="Contact Angle θ [deg]", xlim=(0, max(t_evap_cr, t_evap_ca)*1.05), ylim=(0, 95))
plot!(p3b, sol_cr.t, rad2deg.(sol_cr.theta), label="Constant Radius (CR)", color=:navy, lw=2.5)
plot!(p3b, sol_ca.t, rad2deg.(sol_ca.theta), label="Constant Angle (CA)", color=:darkgreen, lw=2.5, ls=:dash)
plot!(p3b, sol_ss.t, rad2deg.(sol_ss.theta), label="Stick-Slip (SS)", color=:purple, lw=2, ls=:dashdot)
plot!(p3b, sol_cooled.t, rad2deg.(sol_cooled.theta), label="Substrate Cooled (CR)", color=:crimson, lw=2, ls=:dot)

fig3 = plot(p3a, p3b, layout=(2, 1), size=(850, 750))
fig3_path = joinpath(ARTIFACT_DIR, "fig3_sessile_modes_evolution.png")
savefig(fig3, fig3_path)
println("   Saved: $fig3_path")

# FIG 4: Flying vs. Sessile Evaporation Comparison
V0_sessile_mm3 = sol_cr.V[1] * 1e9
p4 = plot(title="Flying Convective Droplet vs. Sessile Droplet Lifetime Comparison",
          xlabel="Time [s]", ylabel="Normalized Volume V(t) / V₀ [-]",
          xlim=(0, max(t_evap_cr, t_evap_cooled)*1.05), ylim=(0, 1.05), legend=:topright)
plot!(p4, t_rig, V_rig ./ V0_mm3, label="Flying Droplet (Rigorous, U=0.5m/s)", color=:navy, lw=3)
plot!(p4, sol_cr.t, sol_cr.V .* 1e9 ./ V0_sessile_mm3, label="Sessile CR (Pinned)", color=:teal, lw=2.5)
plot!(p4, sol_ca.t, sol_ca.V .* 1e9 ./ V0_sessile_mm3, label="Sessile CA (Receding)", color=:darkgreen, lw=2.5, ls=:dash)
plot!(p4, sol_ss.t, sol_ss.V .* 1e9 ./ V0_sessile_mm3, label="Sessile Stick-Slip (SS)", color=:purple, lw=2, ls=:dashdot)
plot!(p4, sol_cooled.t, sol_cooled.V .* 1e9 ./ V0_sessile_mm3, label="Sessile Cooled (Glass Substrate)", color=:crimson, lw=2.5, ls=:dot)

fig4_path = joinpath(ARTIFACT_DIR, "fig4_flying_vs_sessile_lifetime.png")
savefig(p4, fig4_path)
println("   Saved: $fig4_path")

println("\n==================================================")
println("SIMULATION & PLOT GENERATION FINISHED SUCCESSFULLY!")
println("==================================================")
