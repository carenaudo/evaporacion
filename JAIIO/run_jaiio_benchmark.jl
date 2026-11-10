"""
    run_jaiio_benchmark.jl

Benchmark runner for JAIIO 2026 technical report.
Evaluates:
  1. Rigorous Level 3 (Abramzon-Sirignano with parabolic internal temperature profile & Brent root-finding)
  2. Level 2 (Pure convection without Stefan flow blowing correction)
  3. Uniform (Level 3 physics with uniform internal temperature Ts = T_bar, no Brent root-finding)
  4. Isothermal Twb (Iso-Twb, constant temperature wet-bulb droplet)
  5. Modelo 6b Refit (Piecewise Reynolds-scaled empirical surrogate)
  6. Modelo 7 (Multi-variable empirical surrogate with RH and Theta scaling)
  7. Modelo 2b (Optimized d²-law base surrogate)
  8. Modelo 8 Dynamic (Coupled 2-equation volume-temperature dynamic surrogate)

Saves CSV and Markdown summary tables directly into JAIIO/.
"""

using Printf
using Dates
using DataFrames
using CSV
using Statistics
using OrdinaryDiffEq

# Include the Evaporacion package from julia/
const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const JULIA_SRC = joinpath(REPO_ROOT, "julia", "src")
push!(LOAD_PATH, joinpath(REPO_ROOT, "julia"))

using Evaporacion
using Evaporacion.Common
using Evaporacion.Droplet
using Evaporacion.Surrogates

println("="^80)
println("JAIIO 2026: EVALUACIÓN COMPARATIVA DE MODELOS SURROGADOS Y FÍSICOS")
println("Timestamp: ", Dates.now())
println("="^80)

# -----------------------------------------------------------------------------
# 1. Representative Cases Definition
# -----------------------------------------------------------------------------
# Each case tuple: (case_id, d0_m, U_ms, T0_C, T_inf_C, HR)
const CASES = [
    # Baseline test case (Case 1 from benchmark)
    (1, 100e-6, 3.0, 10.0, 50.0, 0.3, "Smallest, fast air, hot, dry"),
    # Standard droplet (Case 13)
    (13, 500e-6, 1.0, 15.0, 30.0, 0.5, "Medium droplet, standard ambient"),
    # Large droplet (Case 22)
    (22, 2000e-6, 1.0, 35.0, 45.0, 0.5, "Large droplet, hot ambient")
]

# -----------------------------------------------------------------------------
# 2. Model Runners
# -----------------------------------------------------------------------------

function run_model_rigorous(d0, U, T0, Tinf, HR; model_variant=:rigorous, t_max=30000.0)
    # Saturation and vapor mass fraction calculation
    twb = stull_twb(Tinf, HR * 100.0)
    psat_inf = psat_water(Tinf + T_ZERO_C)
    pv_inf = HR * psat_inf
    y_inf = surface_mass_fraction_from_psat(pv_inf, P_ATM)

    # Initial droplet mass
    rho_l0 = liquid_density_water(T0 + T_ZERO_C)
    m0 = rho_l0 * (π / 6.0) * (d0^3)

    params = if model_variant === :rigorous
        DropletParams(
            T_inf = Tinf,
            U_rel = U,
            Y_inf = y_inf,
            P_atm = P_ATM,
            profile = :parabolico,
            as1989_blowing = true,
            isothermal = false
        )
    elseif model_variant === :level_2
        DropletParams(
            T_inf = Tinf,
            U_rel = U,
            Y_inf = y_inf,
            P_atm = P_ATM,
            profile = :parabolico,
            as1989_blowing = false, # No Stefan flow blowing correction
            isothermal = false
        )
    elseif model_variant === :uniform
        DropletParams(
            T_inf = Tinf,
            U_rel = U,
            Y_inf = y_inf,
            P_atm = P_ATM,
            profile = :uniforme, # Uniform profile: Ts = T_bar, bypasses Brent
            as1989_blowing = true,
            isothermal = false
        )
    elseif model_variant === :iso_twb
        DropletParams(
            T_inf = Tinf,
            U_rel = U,
            Y_inf = y_inf,
            P_atm = P_ATM,
            profile = :uniforme,
            as1989_blowing = true,
            isothermal = true,
            T_d_fixed = twb # Held fixed at wet-bulb temperature
        )
    else
        error("Unknown model variant: \$model_variant")
    end

    t0 = time_ns()
    sol = solve_droplet(params, [m0, T0], (0.0, t_max); solver=Tsit5(), reltol=1e-6, abstol=1e-9)
    elapsed_s = (time_ns() - t0) / 1e9

    t_evap = sol.t[end]
    n_steps = length(sol.t)
    return (t_evap=t_evap, elapsed_s=elapsed_s, n_steps=n_steps)
end

function run_model_surrogate(d0, U, T0, Tinf, HR; model_type=:modelo_6b_refit, t_max=30000.0)
    d0_mm = d0 * 1e3
    twb = stull_twb(Tinf, HR * 100.0)

    cfg = SurrogateConfig(
        model_type = model_type,
        d0_mm = d0_mm,
        T_inf = Tinf,
        T_wb = twb,
        U = U,
        HR = HR,
        T0 = T0,
        t_end = t_max
    )

    t0 = time_ns()
    sol = solve_surrogate(cfg; solver=Tsit5(), rtol=1e-6, atol=1e-9)
    elapsed_s = (time_ns() - t0) / 1e9

    t_evap = sol.t[end]
    n_steps = length(sol.t)
    return (t_evap=t_evap, elapsed_s=elapsed_s, n_steps=n_steps)
end

# -----------------------------------------------------------------------------
# 3. Main Benchmark Execution
# -----------------------------------------------------------------------------

function run_all_models_single_case(case_tuple)
    cid, d0, U, T0, Tinf, HR, desc = case_tuple
    @printf("\nRunning Case %d: d₀ = %.3f mm, U = %.1f m/s, T₀ = %.1f °C, T_inf = %.1f °C, HR = %d%% (%s)\n",
            cid, d0*1e3, U, T0, Tinf, round(Int, HR*100), desc)
    println("-"^80)

    models_to_test = [
        (:rigorous,         "Rigorous (Level 3)",      :physical,  :rigorous),
        (:level_2,          "Level 2 (Pure Convection)",:physical,  :level_2),
        (:uniform,          "Uniform (No Brent)",      :physical,  :uniform),
        (:iso_twb,          "Iso-Twb (Isothermal)",    :physical,  :iso_twb),
        (:modelo_6b_refit,  "Modelo 6b Refit",         :surrogate, :modelo_6b_refit),
        (:modelo_7,         "Modelo 7 (Multi-Var)",    :surrogate, :modelo_7),
        (:modelo_2b,        "Modelo 2b (d²-law)",      :surrogate, :modelo_2b),
        (:modelo_8_dynamic, "Modelo 8 Dynamic (2-Eq)", :surrogate, :modelo_8_dynamic)
    ]

    results = []

    # Reference lifetime for speedup and error calculation
    ref_time = 0.0
    ref_evap = 0.0

    for (m_id, m_name, m_cat, m_key) in models_to_test
        @printf("  Evaluating %-28s ... ", m_name)

        # Warmup
        if m_cat === :physical
            run_model_rigorous(d0, U, T0, Tinf, HR; model_variant=m_key)
            res = run_model_rigorous(d0, U, T0, Tinf, HR; model_variant=m_key)
        else
            run_model_surrogate(d0, U, T0, Tinf, HR; model_type=m_key)
            res = run_model_surrogate(d0, U, T0, Tinf, HR; model_type=m_key)
        end

        if m_id === :rigorous
            ref_time = res.elapsed_s
            ref_evap = res.t_evap
        end

        speedup = ref_time > 0.0 ? (ref_time / res.elapsed_s) : 1.0
        err_pct = ref_evap > 0.0 ? abs(res.t_evap - ref_evap) / ref_evap * 100.0 : 0.0

        @printf("t_evap: %8.2f s | Time: %7.4f s | Steps: %4d | Speedup: %6.1fx | Err: %5.2f%%\n",
                res.t_evap, res.elapsed_s, res.n_steps, speedup, err_pct)

        push!(results, (
            case_id = cid,
            d0_mm = d0 * 1e3,
            U_ms = U,
            T0_C = T0,
            Tinf_C = Tinf,
            HR = HR,
            model_id = string(m_id),
            model_name = m_name,
            category = string(m_cat),
            t_evap_s = res.t_evap,
            elapsed_s = res.elapsed_s,
            n_steps = res.n_steps,
            speedup = speedup,
            error_pct = err_pct
        ))
    end

    return DataFrame(results)
end

# Run the benchmark for each case
all_dfs = DataFrame[]
for c in CASES
    df_case = run_all_models_single_case(c)
    push!(all_dfs, df_case)
end

final_df = vcat(all_dfs...)

# Save results CSV
out_csv = joinpath(@__DIR__, "jaiio_benchmark_results.csv")
CSV.write(out_csv, final_df)
@printf("\nSaved benchmark results to: %s\n", out_csv)

# Print Summary Table
println("\n" * "="^80)
println("JAIIO 2026 BENCHMARK EXECUTION SUMMARY")
println("="^80)
println(final_df[:, [:case_id, :model_id, :t_evap_s, :elapsed_s, :n_steps, :speedup, :error_pct]])

println("\nDone! All JAIIO models calculated successfully.")
