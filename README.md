# Evaporacion.jl – Droplet Evaporation Simulation Suite in High-Performance Julia

`Evaporacion.jl` is a high-performance, physics-grounded simulation library written natively in **Julia** for simulating water droplet evaporation, transient cooling, and drying kinetics across both **in-flight** (falling / suspended) and **sessile** (deposited on substrates) configurations.

---

## 1. Core Modules & Physical Capabilities

- **`Evaporacion.Droplet`**: In-flight spherical droplets with coupled mass and energy ODEs, implicit surface temperature solving via Brent's method, Abramzon & Sirignano (1989) blowing closures, non-unitary Lewis numbers, dynamic free-fall terminal velocity, and Kronig-Brink internal liquid circulation. Supports `:d2_classic`, `:d2_stefan`, `:d2_convective`, and `:abramzon_sirignano` variants.
- **`Evaporacion.Sessile`**: Sessile droplets sitting on impermeable substrates across all four canonical modes: Constant Contact Radius (`:CR`), Constant Contact Angle (`:CA`), Stick-Slide (`:SS`), and Stick-Jump (`:SJ`). Incorporates exact Mehler-Fock toroidal integral evaluation of $g(\theta)$, Popov contact-line singular flux, substrate conductive cooling, internal lubrication flow velocity fields $(\bar{u}, u, w)$, and coffee-ring deposit scaling (Wilson & D'Ambrosio 2023).
- **`Evaporacion.RadialPDE`**: 1D moving-boundary radial conduction PDE solved using the Method of Lines on an immobilized Landau coordinate $\xi = r/R(t)$, featuring Marangoni interfacial stresses, Rayleigh natural convection, and Langmuir surfactant adsorption/desorption.
- **`Evaporacion.Surrogates`**: Calibrated Reynolds-segmented 1-equation and 2-equation volume and temperature surrogates (`modelo_2b` to `modelo_8_dynamic`), providing up to 600× speedups over the full non-linear ODE solver.
- **`Evaporacion.Neural`**: SciML Neural ODE surrogates in Lux.jl, including pure MLP Neural ODEs (`Option3bNeuralODE`) and hybrid physical + regime-attention transformer ODEs (`Option3cHybridODE`).

---

## 2. Directory Structure

```
evaporacion/
├── julia/               # Core Julia package
│   ├── Project.toml     # Package dependencies (OrdinaryDiffEq, Lux, QuadGK, Plots, etc.)
│   ├── Manifest.toml    # Exact dependency lockfile
│   ├── src/             # Package source code (Common, Droplet, Sessile, RadialPDE, Surrogates, Neural)
│   ├── test/            # Test suite (runtests.jl)
│   └── evaporacion_gota.jl # Package loading entry point
├── papers/              # 60 Curated literature PDFs and standardized mathematical summaries
│   ├── *.pdf            # Cleanly renamed author_year.pdf files
│   └── *_summary.md     # Standardized summaries (equations, nomenclature, validity limits)
├── clicap/              # CLICAP 2025 paper reproduction suite & forensic reports
│   ├── reproduce_paper_cases.jl  # In-flight & sessile reproduction script
│   ├── generate_figure2.jl       # Sessile evaporation modes comparison script
│   ├── run_case_study.jl         # Single-case in-flight vs sessile comparison
│   └── *.png, *.csv, *.md        # Reproduction figures, tables, and reports
├── JAIIO/               # JAIIO 2026 surrogate benchmark reports & validation studies
│   ├── REPORT_JAIIO.md           # Technical report on surrogate accuracy and speedups
│   ├── validation_report.md      # Detailed 2520-case validation report
│   └── validation_report_regimes.md # Regime-separated (diffusion vs convection) validation
├── docs/                # Architectural guides & technical documentation
│   ├── KNOWLEDGE_GRAPH.md        # Low-token call graph & symbol-to-line index
│   ├── MODELS_COMPARISON_AND_CHECKLIST.md # Model comparison & missing-features checklist
│   ├── PAPERS_CATALOG.md         # Master index of the 60 literature publications
│   ├── RIGOROUS_MODEL_SPEC.md    # Level 3 reference model specification
│   ├── AS1989_EQUATIONS_AND_WILKE.md # Abramzon-Sirignano equations & Wilke mixing rules
│   └── BT_NEGATIVO_DIAGNOSTICO.md # Diagnostic of surface energy balance & negative Spalding BT
└── legacy/              # Archived legacy Python code, scratch scripts, and outputs (to be removed)
```

---

## 3. Quick Start (Julia)

### 3.1 Activate Environment & Run Tests

```bash
julia --project=julia -e "using Pkg; Pkg.instantiate(); Pkg.test()"
```

### 3.2 Run a Droplet Evaporation Simulation

```julia
using Evaporacion

# Configure parameters for an evaporating droplet
params = DropletParams(
    T_inf = 35.0,    # Ambient temperature [°C]
    U_rel = 0.5,     # Relative air speed [m/s]
    Y_inf = 0.015,   # Ambient vapor mass fraction
    profile = :parabolico
)

# Initial conditions: mass for a 1 mm droplet, initial temperature 22 °C
m0 = liquid_density_water(22.0 + 273.15) * (π / 6.0) * (0.001^3)
u0 = [m0, 22.0]
tspan = (0.0, 400.0)

# Solve ODE
sol = solve_droplet(params, u0, tspan)
df = postprocess_trajectory(sol, params)
println("Evaporation lifetime: ", round(df.t[end], digits=2), " s")
```

### 3.3 Run Sessile Droplet Lifetimes

```julia
using Evaporacion

# Configure sessile droplet in Constant Radius (pinned) mode
cfg = SessileDropletConfig(
    R0 = 1.0e-3,           # Initial contact radius [m]
    theta0 = 0.5 * π,      # Contact angle [rad] (90°)
    mode = :CR
)

t_cr = lifetime_cr(cfg)
sol_cr = solve_sessile(cfg)
println("Sessile CR lifetime: ", round(t_cr, digits=2), " s")
```

### 3.4 Run CLICAP Paper Cases

```bash
julia --project=julia clicap/reproduce_paper_cases.jl
```

---

## 4. Documentation Index

- **[Knowledge Graph (`docs/KNOWLEDGE_GRAPH.md`)](file:///d:/programacion/evaporacion/docs/KNOWLEDGE_GRAPH.md)**: Fast, token-efficient architecture map and line-level code locator.
- **[Models Comparison & Checklist (`docs/MODELS_COMPARISON_AND_CHECKLIST.md`)](file:///d:/programacion/evaporacion/docs/MODELS_COMPARISON_AND_CHECKLIST.md)**: Detailed feature-by-feature matrix comparing `Evaporacion.jl` against the 60 papers in `papers/`.
- **[Master Papers Catalog (`docs/PAPERS_CATALOG.md`)](file:///d:/programacion/evaporacion/docs/PAPERS_CATALOG.md)**: Index and classification of all 60 publications.
- **[JAIIO Surrogate Studies (`JAIIO/`)](file:///d:/programacion/evaporacion/JAIIO/)**: Performance benchmarks, Pareto charts, and speedup evaluations.
