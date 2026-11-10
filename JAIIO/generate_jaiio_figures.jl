"""
    generate_jaiio_figures.jl

Generates benchmark comparison figures for the JAIIO 2026 report using Plots.jl:
  1. Execution time barplot / scatter across models
  2. Speedup comparison relative to rigorous Level 3
  3. Pareto trade-off curve (Execution time vs MAPE/Error %)

Saves figures directly into JAIIO/ as PNG files.
"""

using CSV
using DataFrames
using Plots
using Printf
using Statistics

# Configure GR headless environment for clean plot generation
ENV["GKSwstype"] = "100"

const CSV_PATH = joinpath(@__DIR__, "jaiio_benchmark_results.csv")
if !isfile(CSV_PATH)
    error("Benchmark CSV not found at: \$CSV_PATH. Run run_jaiio_benchmark.jl first.")
end

df = CSV.read(CSV_PATH, DataFrame)

# Group by model to compute average speedup, elapsed time, and error %
gdf = combine(groupby(df, :model_name),
    :elapsed_s => mean => :mean_elapsed_s,
    :speedup => mean => :mean_speedup,
    :error_pct => mean => :mean_error_pct,
    :n_steps => mean => :mean_steps
)

# Sort models by execution speed (slowest to fastest)
sort!(gdf, :mean_elapsed_s, rev=true)

println("Aggregated JAIIO Model Metrics:")
println(gdf)

# -----------------------------------------------------------------------------
# 1. Figure 1: Mean Execution Time (Log Scale)
# -----------------------------------------------------------------------------
p1 = bar(
    gdf.model_name,
    gdf.mean_elapsed_s,
    yscale = :log10,
    ylabel = "Tiempo de Ejecución Medio [s] (log₁₀)",
    title = "JAIIO 2026: Tiempo Computacional por Modelo",
    legend = false,
    color = :dodgerblue,
    xrotation = 35,
    margin = 5Plots.mm,
    size = (800, 500)
)
savefig(p1, joinpath(@__DIR__, "jaiio_exec_time_barplot.png"))
println("Saved: jaiio_exec_time_barplot.png")

# -----------------------------------------------------------------------------
# 2. Figure 2: Mean Speedup relative to Rigorous Level 3
# -----------------------------------------------------------------------------
# Filter out the rigorous model itself for the speedup barplot
gdf_surr = filter(row -> row.model_name != "Rigorous (Level 3)", gdf)
sort!(gdf_surr, :mean_speedup)

p2 = bar(
    gdf_surr.model_name,
    gdf_surr.mean_speedup,
    ylabel = "Factor de Aceleración (Speedup) vs Riguroso",
    title = "JAIIO 2026: Speedup Relativo al Modelo Riguroso (Nivel 3)",
    legend = false,
    color = :forestgreen,
    xrotation = 35,
    margin = 5Plots.mm,
    size = (800, 500)
)
savefig(p2, joinpath(@__DIR__, "jaiio_speedup_barplot.png"))
println("Saved: jaiio_speedup_barplot.png")

# -----------------------------------------------------------------------------
# 3. Figure 3: Pareto Front (Execution Time vs Error %)
# -----------------------------------------------------------------------------
p3 = scatter(
    gdf.mean_elapsed_s,
    gdf.mean_error_pct,
    xscale = :log10,
    xlabel = "Tiempo de Ejecución Medio [s] (escala logarítmica)",
    ylabel = "Error Medio Respecto a Riguroso [%]",
    title = "JAIIO 2026: Frontera de Pareto (Precisión vs Velocidad)",
    markersize = 8,
    color = :crimson,
    label = "Modelos Evaluados",
    margin = 5Plots.mm,
    size = (800, 500)
)

for row in eachrow(gdf)
    annotate!(row.mean_elapsed_s * 1.15, row.mean_error_pct, text(row.model_name, 8, :left))
end

savefig(p3, joinpath(@__DIR__, "jaiio_pareto_accuracy_vs_speed.png"))
println("Saved: jaiio_pareto_accuracy_vs_speed.png")

println("\nAll JAIIO figures generated successfully!")
