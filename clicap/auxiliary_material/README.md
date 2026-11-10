# Supplementary Material: Multi-Model Droplet Evaporation Benchmark & Reproduction Guide

This directory constitutes the official **Supplementary / Auxiliary Material** accompanying the publication on droplet evaporation physics and multi-model benchmarking.

All computational models are implemented in native, high-performance Julia via **`Evaporacion.jl`**.

---

## 1. Directory Structure

```text
clicap/auxiliary_material/
├── README.md                                  # This execution and guide document
├── CLICAP_MODELS_COMPARISON_REPORT.md         # Comprehensive scientific benchmark report
├── scripts/
│   ├── run_all_cases.jl                       # Master script: runs all cases, generates all CSVs & PNGs
│   ├── run_inflight_cases.jl                  # In-flight droplet simulations across 9 cases, 4 diameters & 10 models
│   ├── run_sessile_cases.jl                   # Sessile droplet simulations across 9 cases, 3 angles & 8 models
│   ├── compute_surface_areas.jl               # Computes time-averaged mean and median surface areas
│   └── run_single_case.jl                     # Interactive / CLI runner for arbitrary custom scenarios
├── data/
│   ├── table1_environmental_cases.csv         # Matrix of 9 environmental cases (T, HR, Twb, psat, k, d_min)
│   ├── table2_inflight_cases1_and_7.csv       # Comparison with paper Table 2 (Cases 1 & 7)
│   ├── table3_sessile_paper_comparison.csv    # Comparison with paper Table 3 (Sessile modes)
│   ├── comparison_inflight_all_cases.csv      # Complete in-flight dataset (360 simulation runs)
│   ├── comparison_sessile_all_cases.csv       # Complete sessile dataset (27 parameter combinations)
│   ├── comparison_cross_domain_flying_vs_sessile.csv # Direct cross-domain flying vs sessile comparison
│   ├── table_mean_median_surface_areas_summary.csv   # Summary table of mean and median surface areas
│   ├── surface_areas_inflight_all_cases.csv   # Full surface areas dataset for flying droplets (360 rows)
│   └── surface_areas_sessile_all_cases.csv    # Full surface areas dataset for sessile droplets (162 rows)
└── figures/
    ├── fig_inflight_lifetimes_by_diameter.png # Evaporation lifetime vs d0 for Cases 1 & 7
    ├── fig_inflight_K_initial_vs_HR.png       # Initial evaporation rate K(t=0) vs RH (Extended Fig 3)
    ├── fig_inflight_extinction_dmin_vs_weather.png # Extinction cut-off diameter d_min vs RH and T
    ├── fig_sessile_lifetimes_vs_HR.png        # Sessile modes lifetime vs RH (Extended Fig 4)
    ├── fig_sessile_modes_vs_contact_angle.png # Sessile modes lifetime vs contact angle theta0
    └── fig_flying_vs_sessile_all_cases.png    # Global cross-domain lifetime across all 9 cases
```

---

## 2. Requirements & Environment Setup

All scripts require **Julia 1.10+** (or Julia 1.11+).

To activate the project environment with all required dependencies:

```bash
julia --project=julia -e "using Pkg; Pkg.instantiate()"
```

Dependencies installed in `julia/Project.toml`:
- `OrdinaryDiffEq.jl` (High-order adaptive ODE solvers, Tsit5 / Rodas5P)
- `DataFrames.jl` & `CSV.jl` (Tabular data generation and export)
- `Plots.jl` with `GR` backend (Publication-quality figures)
- `QuadGK.jl` (Adaptive Gauss-Kronrod numerical quadrature)
- `ArgParse.jl` (Command-line interface argument parsing)

---

## 3. How to Reproduce All Results

To run the complete benchmark suite from scratch (generating all 6 CSV tables and 6 figures):

```bash
julia --project=julia clicap/auxiliary_material/scripts/run_all_cases.jl
```

Execution time: ~50 seconds for 360 non-linear ODE integrations, 27 analytical quadrature integrals, and high-resolution figure rendering.

---

## 4. How to Run Individual Model & Scenario Scripts

### A. Run All In-Flight Droplet Cases
Simulates all 9 environmental cases for $d_0 \in \{50, 200, 600, 1200\}\ \mu\text{m}$ across all 10 in-flight models:
```bash
julia --project=julia clicap/auxiliary_material/scripts/run_inflight_cases.jl
```

### B. Run All Sessile Droplet Cases
Simulates all 9 environmental cases for $d_0^* = 200\ \mu\text{m}$ across contact angles $\theta_0 \in \{20^\circ, 40^\circ, 80^\circ\}$ and all 8 sessile models:
```bash
julia --project=julia clicap/auxiliary_material/scripts/run_sessile_cases.jl
```

### C. Run Any Custom Single Case via CLI (`run_single_case.jl`)

The `run_single_case.jl` script provides a flexible CLI for running any specific scenario or model:

#### In-Flight Droplet Examples:
```bash
# Run Abramzon-Sirignano for Case 1 (10 °C, 20% HR) with d0 = 200 µm:
julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain inflight --case 1 --d0 200 --model abramzon_sirignano

# Run Classical D²-law for Case 7 (40 °C, 20% HR) with d0 = 600 µm:
julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain inflight --case 7 --d0 600 --model d2_classic

# Run Zhou et al. (2017) 3rd-order profile model for Case 4 (25 °C, 20% HR):
julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain inflight --case 4 --d0 200 --model zhou_cubic

# Run Vlasov (2021) Knudsen-layer kinetic model for a 50 µm microdroplet:
julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain inflight --case 1 --d0 50 --model vlasov_kinetic

# Run with custom weather parameters (e.g. T = 32.5 °C, HR = 45%):
julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain inflight --T 32.5 --HR 45.0 --d0 250 --model abramzon_sirignano
```

#### Sessile Droplet Examples:
```bash
# Run Constant Contact Radius (CR / Pinned mode) for Case 7 with theta0 = 40°:
julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain sessile --case 7 --theta0 40 --model CR

# Run Constant Contact Angle (CA / Receding mode) for Case 1 with theta0 = 80°:
julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain sessile --case 1 --theta0 80 --model CA

# Run Stick-Slide (SS mode) for Case 4 with theta0 = 20°:
julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain sessile --case 4 --theta0 20 --model SS

# Run Substrate Conductive Cooling on Glass:
julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain sessile --case 7 --theta0 40 --model substrate_cooling

# Run Cooperative Multi-Droplet Vapor Shielding:
julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain sessile --case 5 --theta0 40 --model vapor_shielding
```

#### CLI Available Options Table:
| Option | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--domain` | String | `"inflight"` | Evaporation domain: `"inflight"` or `"sessile"` |
| `--case` | Int | `0` | Standard CLICAP case (1 to 9). Overrides `--T` and `--HR` |
| `--T` | Float64 | `25.0` | Ambient temperature in Celsius (used when `--case 0`) |
| `--HR` | Float64 | `50.0` | Ambient relative humidity in % (used when `--case 0`) |
| `--P_atm` | Float64 | `101325.0` | Atmospheric pressure in Pa |
| `--d0` | Float64 | `200.0` | Initial diameter in $\mu\text{m}$ (or equivalent volume for sessile) |
| `--theta0`| Float64 | `40.0` | Initial contact angle in degrees (sessile droplets) |
| `--model` | String | `"abramzon_sirignano"` | Model key (see list below) |
| `--U_rel` | Float64 | `0.0` | Relative air velocity in m/s (0 = auto terminal velocity) |
| `--output`| String | `""` | Optional CSV destination file |

Supported model identifiers:
- **In-Flight**: `d2_classic`, `d2_stefan`, `d2_convective`, `d2_revised`, `abramzon_sirignano`, `abramzon_uniform`, `zhou_cubic`, `vlasov_kinetic`, `surrogate_6b`, `surrogate_8dyn`
- **Sessile**: `CA`, `CR`, `SS`, `SJ`, `thin_cr`, `thin_ca`, `picknett`, `substrate_cooling`, `vapor_shielding`

---

## 5. Mean & Median Surface Area in Time Computation

To calculate the time-averaged mean surface area $\bar{S} = \frac{1}{t_{\text{evap}}} \int_0^{t_{\text{evap}}} S(t)\mathrm{d}t$ and temporal median surface area $S_{\text{med}} = S(t = t_{\text{evap}}/2)$ across all 9 cases and models:

```bash
julia --project=julia clicap/auxiliary_material/scripts/compute_surface_areas.jl
```

This exports:
- `data/table_mean_median_surface_areas_summary.csv`: Master publication table comparing $S_0$, mean, and median surface areas between in-flight and sessile modes for all 9 environmental cases.
- `data/surface_areas_inflight_all_cases.csv`: Detailed dataset for 360 in-flight parameter combinations ($d_0 \in \{50, 200, 600, 1200\}\ \mu\text{m}$, 10 models).
- `data/surface_areas_sessile_all_cases.csv`: Detailed dataset for 162 sessile parameter combinations ($\theta_0 \in \{20^\circ, 40^\circ, 80^\circ\}$, 6 modes).

---

## 6. Dynamic Volume, V/S Ratio, & Computational Performance Benchmark

To calculate the time-averaged mean volume $\bar{V}$, median volume $V_{\text{med}}$, the dynamic hydraulic thickness ratio $(V/S)$, and perform high-resolution timing analysis (1 warm-up discarded, 5 timed runs per model):

```bash
julia --project=julia clicap/auxiliary_material/scripts/compute_volumes_and_timings.jl
```

This exports:
- `data/table_mean_median_volume_surface_summary.csv`: Cross-domain comparison of $\bar{V}$, $\bar{S}$, $\bar{V}/\bar{S}$, and $(V/S)$ ratios for in-flight vs sessile droplets ($d_0 = 200\ \mu\text{m}$, $\theta_0 = 40^\circ$) across all 9 environmental cases.
- `data/volumes_surface_ratios_inflight_all_cases.csv`: Detailed dataset for all 360 in-flight configurations.
- `data/volumes_surface_ratios_sessile_all_cases.csv`: Detailed dataset for all 162 sessile configurations.
- `data/table_model_computation_times_summary.csv`: Average calculation times [ms] across all 9 cases and relative speedup factors vs Abramzon-Sirignano.
- `data/calculation_times_all_cases.csv`: Detailed per-case execution times (mean, median, min, max, std).
