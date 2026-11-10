# Walkthrough: Port of Evaporation Suite to Julia (`Evaporacion.jl`)

The complete Python evaporation suite has been ported, audited against primary literature, and modernized as a high-performance Julia package named **`Evaporacion.jl`**.

All legacy polynomial approximations and compatibility shims have been retired, standardizing strictly on canonical **High-Fidelity (Level 3)** physical property correlations and literature closures.

---

## 1. Package Architecture

```
julia/
├── Project.toml                      # Modern Pkg specification with full dependencies & compats
├── evaporacion_gota.jl               # Canonical top-level runner script delegating to Evaporacion.CLI
├── src/
│   ├── Evaporacion.jl                # Top-level umbrella module and unified exports
│   │
│   ├── common/                       # Component 1: Shared Core
│   │   ├── Common.jl
│   │   ├── constants.jl              # Universal constants (R, M_a, M_w, σ_SB, P_atm)
│   │   ├── properties.jl             # Buck, Kell 5th-order, Watson, Sutherland, Fuller, Stull
│   │   ├── correlations.jl           # Frössling, Clift, A&S blowing closures, Kronig-Brink
│   │   └── geometry.jl               # Spherical & Sessile geometry, exact g(θ) & O(1) spline
│   │
│   ├── droplet/                      # Component 2: Rigorous Spherical Model
│   │   ├── Droplet.jl
│   │   ├── types.jl                  # DropletParams & DropletAlgebraic (37 diagnostic fields)
│   │   ├── algebraic.jl              # Dual-compatible Brent root-finding for surface Ts
│   │   ├── solver.jl                 # OrdinaryDiffEq ODEs (Rodas5P, RadauIIA5) + mass callback
│   │   └── postprocess.jl            # Vectorized diagnostic trajectory DataFrame builder
│   │
│   ├── surrogates/                   # Component 3: Calibrated Empirical & Coupled Models
│   │   ├── Surrogates.jl
│   │   ├── algebraic.jl              # Modelos 2b, 4, 5, 6, 6b, 6b_refit, 7, 8, 8_dynamic
│   │   ├── solver.jl                 # 1-variable (V) and 2-variable (V, T) surrogate solvers
│   │   └── ensemble.jl               # Multithreaded parallel sweep runner via EnsembleProblem
│   │
│   ├── sessile/                      # Component 4: Sessile Droplet Suite (Wilson 2023)
│   │   ├── Sessile.jl
│   │   ├── flux.jl                   # Popov (2005) global/local flux, substrate cooling, Deegan
│   │   ├── solver.jl                 # CR, CA, SS, and SJ modes + exact analytical lifetimes
│   │   └── internal_deposition.jl   # Eqs. (24)–(26) velocity fields & coffee-ring accumulation
│   │
│   ├── radial_pde/                   # Component 5: 1D Moving-Boundary Radial PDE
│   │   ├── RadialPDE.jl
│   │   └── solver.jl                 # Immobilized Landau coordinate ξ ∈ [0, 1] MOL solver
│   │
│   ├── neural/                       # Component 6: SciML Neural ODE Surrogates
│   │   ├── Neural.jl
│   │   ├── option3b_mlp.jl           # Option 3b Pure MLP Neural ODE in Lux.jl
│   │   ├── option3c_hybrid_attn.jl   # Option 3c Hybrid Regime-Attention Neural ODE
│   │   ├── checkpoint_io.jl          # Zero-dependency PyTorch NPZ checkpoint importer
│   │   └── train.jl                  # SciMLSensitivity & Optimisers adjoint training loop
│   │
│   ├── plotting.jl                   # Diagnostic trajectory & profile visualizers (Plots.jl)
│   └── cli.jl                        # ArgParse CLI supporting all models & output formats
│
├── test/
│   └── runtests.jl                   # Full regression and physics parity test suite
└── archive/
    └── legacy_julia_scratch/         # Archived historical scratch scripts
```

---

## 2. Mathematical Alignments & Bug Fixes

### A. Abramzon & Sirignano (1989) Hill-Vortex Internal Circulation
- **Paper Reference**: Section 4, p. 1610, Eqs. (30) & (32).
- **Reynolds Exponent Correction**:
  $$C_F = \frac{12.69}{Re^{2/3} (1 + B_M)}, \qquad U_s = \frac{1}{32} U \, Re \left(\frac{\mu_g}{\mu_l}\right) C_F = \frac{\rho_g U^2 r_S}{16 \mu_l} C_F$$
  *(Corrected from older $Re^{-1/3}$ code mismatches to match the primary literature exactly).*

### B. Wilson & D'Ambrosio (2023) Sessile Cap Evolution
- **Exact Analytical Identity**:
  $$V(R, \theta) = \frac{\pi R^3}{3} \frac{\sin\theta (2+\cos\theta)}{(1+\cos\theta)^2} \implies \frac{\partial V}{\partial \theta} = \frac{\pi R^3}{(1+\cos\theta)^2}$$
- **Contact Line Pinned Mode (CR Mode)**:
  $$\frac{d\theta}{dt} = - \frac{D \Delta c_{\text{eff}}}{\rho_l R_0^2} g(\theta)$$
  The $(1+\cos\theta)^2$ term cancels identically on both sides, eliminating catastrophic cancellation.
- **Analytical Lifetimes**:
  Implemented exact closed-form and integral formulas from Tables 1 & 2 for all four evaporation regimes:
  - $t_{\text{CR}} = \frac{\rho_l R_0^2}{D \Delta c} \int_0^{\theta_0} \frac{1}{g(\theta)} d\theta$
  - $t_{\text{CA}} = \frac{\rho_l R_0^2}{2 D \Delta c} \frac{\sin\theta_0 (2+\cos\theta_0)}{g(\theta_0)}$
  - $t_{\text{SS}} = \frac{\rho_l R_0^2}{2 D \Delta c} \left( \int_{\theta^*}^{\theta_0} \frac{2}{g(\theta)} d\theta + \frac{\sin\theta^* (2+\cos\theta^*)}{g(\theta^*)} \right)$
  - $t_{\text{SJ}} = \frac{\rho_l R_0^2}{D \Delta c} \left( \int_{\theta_{\min}}^{\theta_{\max}} \frac{1}{g(\theta)} d\theta \right) \frac{1}{1 - \Lambda^{2/3}}$

### C. 1D Radial Moving-Boundary PDE Corrections
- **Convective Coefficient**: Replaced the hardcoded constant $h_{\text{conv}} = 10.0$ with the rigorous dynamically computed Nusselt correlation $h_{\text{conv}} = Nu \frac{k_g}{2 R}$.
- **Temperature Units**: Eliminated mixed Celsius/Kelvin additions in film temperature calculations.
- **Raoult / Solute Mass Fraction**: Converted solute mole fractions to proper mass fractions before evaluating $B_M$.

---

## 3. Verification & Parity Results

Standard-library ground truth reference calculations verified the following canonical physical properties:
- **Saturation Vapor Pressure $P_{\text{sat}}(25^\circ\text{C})$**: $3168.531412 \text{ Pa}$ (Buck 1981)
- **Liquid Density $\rho_l(25^\circ\text{C})$**: $997.047958 \text{ kg/m}^3$ (Kell 5th-order)
- **Latent Heat $\Delta H_{\text{ev}}(25^\circ\text{C})$**: $2.474383\times 10^6 \text{ J/kg}$ (Watson 1943)
- **Latent Heat at Normal Boiling $\Delta H_{\text{ev}}(100^\circ\text{C})$**: $2.257000\times 10^6 \text{ J/kg}$
- **Air Viscosity $\mu_g(25^\circ\text{C})$**: $1.837149\times 10^{-5} \text{ Pa}\cdot\text{s}$ (Sutherland)
- **Water Vapor Diffusivity $D_{ab}(25^\circ\text{C})$**: $2.508204\times 10^{-5} \text{ m}^2/\text{s}$ (Fuller et al. 1966)
- **Wet Bulb Temperature $T_{\text{wb}}(35^\circ\text{C}, 40\% \text{ RH})$**: $24.514154^\circ\text{C}$ (Stull 2011)

---

## 4. Usage Instructions

### Running via Top-Level CLI:
```bash
# Rigorous spherical droplet
julia --project=julia julia/evaporacion_gota.jl --model droplet --d0 1.0 --Tinf 35.0 --HR 40.0 --U 0.5 --profile parabolico

# Calibrated surrogate (e.g. Modelo 6b Refit or Modelo 8 Dynamic)
julia --project=julia julia/evaporacion_gota.jl --model surrogate --surrogate-type modelo_6b_refit --d0 1.0 --Tinf 35.0 --HR 40.0

# Sessile droplet (e.g. Stick-Slip mode)
julia --project=julia julia/evaporacion_gota.jl --model sessile --sessile-mode SS --d0 1.0 --t-end 1500.0

# 1D moving-boundary radial PDE
julia --project=julia julia/evaporacion_gota.jl --model radial_pde --d0 1.0 --Tinf 35.0 --HR 40.0
```

### Running in Julia Code:
```julia
using Evaporacion

# 1. Rigorous Droplet
params = DropletParams(T_inf=35.0, U_rel=1.0, Y_inf=0.01, profile=:parabolico)
m0 = liquid_density_water(295.15) * (π / 6.0) * (0.001^3)
sol = solve_droplet(params, [m0, 22.0], (0.0, 100.0))
df = postprocess_trajectory(sol, params)

# 2. Sessile Droplet
cfg_sessile = SessileDropletConfig(R0=1e-3, theta0=0.5*π, mode=:CR)
sol_sessile = solve_sessile(cfg_sessile)
println("Analytical evaporation time: ", lifetime_cr(cfg_sessile), " s")

# 3. Parallel Ensemble Parametric Sweep
configs = [SurrogateConfig(d0_mm=d, T_inf=T) for d in [0.5, 1.0, 2.0] for T in [25.0, 35.0, 45.0]]
sweep_df = run_parametric_sweep(configs)
```
