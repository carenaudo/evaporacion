# Evaporation Models Comparison & Missing Features Checklist

This document provides a comparative analysis between the models implemented in **`Evaporacion.jl`** and the physical/mathematical models formulated across the **60 publications** in [`papers/`](file:///d:/programacion/evaporacion/papers).

---

## 1. Inventory of Evaporation Models in `Evaporacion.jl`

| Julia Submodule | Model Identifier | Governing Differential Equations | Key Closures & Physics | Primary Reference in `papers/` |
| :--- | :--- | :--- | :--- | :--- |
| **`Droplet`** | `:d2_classic` | $\frac{dm}{dt} = -4\pi r_s \rho_g D_{ab} (Y_s - Y_\infty)$ | Pure diffusive $D^2$-law ($\text{Sh}=2$), optional fixed $T_d = T_{\text{wb}}$ (Stull) | Langmuir (1918), Godsave (1953), Renaudo et al. (2025) Model a |
| **`Droplet`** | `:d2_stefan` | $\frac{dm}{dt} = -4\pi r_s \rho_g D_{ab} \ln(1 + B_M)$ | Diffusive $D^2$-law with radial Stefan flow ($B_M = \frac{Y_s - Y_\infty}{1 - Y_s}$) | Spalding (1953), Akkuş (2020), Renaudo et al. (2025) Model b |
| **`Droplet`** | `:d2_convective` | $\frac{dm}{dt} = -2\pi r_s \rho_g D_{ab} \text{Sh}_0 (Y_s - Y_\infty)$ | Frössling convection ($\text{Sh}_0 = 2 + 0.553 \text{Re}^{1/2}\text{Sc}^{1/3}$), dynamic $v_{\text{term}}(d)$ | Kinzer & Gunn (1951), Holterman (2003), Renaudo et al. (2025) Model c |
| **`Droplet`** | `:abramzon_sirignano` (Level 3) | $u = [m, \bar{T}_l]$; root-find for $T_s$: $q_g - \dot{m}\Delta H_{\text{ev}} - q_l = 0$ | Coupled $B_M, B_T = (1+B_M)^\phi - 1$, blowing factors $F(B)$, non-unitary Lewis, Kronig-Brink circulation $\chi_l(\text{Pe}_l)$ | Abramzon & Sirignano (1989), Sirignano (2010), Sazhin (2006) |
| **`Sessile`** | Modes `:CR`, `:CA`, `:SS`, `:SJ` | Exact lifetimes + ODEs for $\theta(t)$ or $R(t)$ | Toroidal Mehler-Fock integral $g(\theta)$ (exact & fast table), Popov/Deegan local flux, Picknett-Bexon polynomials, shallow thin-droplet limits, substrate cooling $\Delta c_{\text{eff}}$, multi-droplet vapor shielding | Picknett & Bexon (1977), Stauber et al. (2015), Wilson & D'Ambrosio (2023), Gelderblom (2022), Pan & Wang (2022) |
| **`Kinetic`** | `:vlasov2021`, `:zhao2023`, `:akkus2020` | Mixed kinetic-diffusion ODE $\frac{dm}{dt}$ with Knudsen transition | Hertz-Knudsen-Schrage jump $I(Kn, \beta)$, kinetic length $l_k = \frac{4D}{\alpha \bar{v}}$, analytical $R(t)$, lifetime extension $\tau/\tau_{\text{Maxwell}} = 1 + 2l_k/R_0$, kinetic Péclet $\text{Pe}_k$, Kelvin curvature $\exp(\frac{2\sigma M_v}{\rho_l R_g T r})$, Stefan blowing $\Theta_{\text{Stefan}}$ | Vlasov (2021), Zhao & Nadal (2023), Akkuş (2020), Sazhin (2006) |
| **`Multicomponent`** | `solve_stefan_maxwell_exact`, `solve_multicomponent_droplet` | Coupled multicomponent ODEs $[m_i(t), T_d(t)]$ | Exact matrix exponential $e^{\mathbf{A}}(\mathbf{A}\mathbf{W}_\infty + \mathbf{B}) - (\mathbf{A}\mathbf{W}_s + \mathbf{B}) = \mathbf{0}$, Ordoubadi Stefan-Fuchs convective fluxes $\dot{m}_i, \dot{m}_t$, Raoult vapor pressure depression | Tonini & Cossali (2015), Ordoubadi et al. (2019) |
| **`Saline`** | `solve_saline_droplet` | Saline/colloidal evaporation ODEs $[m_w(t), m_s, T_d(t)]$ | Raoult-van 't Hoff vapor depression $p_{vs} = p_{\text{sat}}(1 - i Z_s)$, time-varying evaporation constant $\beta(t)$, crystallization time $t_p$, Vehring solute Péclet $\text{Pe}_j = \frac{\kappa}{8 D_j}$ & surface enrichment $E_j$, Gerken fractal agglomerate packing limit $\phi_{\text{max}} = (R_g/r_{np})^{d_f-3}\frac{\pi}{6}$, crust locking diameter $d_{\text{lock}} = d_0 \phi_{s,0}^{1/3}$ | Maruyama & Hasegawa (2020), Gerken et al. (2014), Vehring et al. (2007) |
| **`Common`** | In-flight extinction & spray regimes | Algebraic closures & regime classifications | In-flight evaporation parameter $k = \frac{1.25\times 10^9}{T_c - T_{bh}}$, cut-off extinction diameter $d_{\text{min}} = [\frac{54\mu_g H}{k(\rho_d-\rho_g)g}]^{1/4}$, landing diameter $d_{\text{dep}}(d_0)$, Chiu group evaporation number $G$ | Renaudo et al. (2022, 2023), Chiu (Sirignano 2010), Xu & Xu (2025) |
| **`RadialPDE`** | `RadialPDEConfig` | $y = [R(t), T(\xi_1..\xi_N), m_w, \Gamma]$ | 1D moving-boundary conduction on Landau grid $\xi = r/R(t)$, natural convection $\text{Sh}(\text{Ra})$, Marangoni flow, surfactants | Sazhin (2006), Niazmand et al. (1994), Maruyama (2020), Gelderblom et al. (2022) |
| **`Surrogates`** | 1-Eq: `modelo_2b..7`<br/>2-Eq: `modelo_8(_dynamic)` | $\frac{dV}{dt} = -C \sqrt{S}(T_\infty - T_{\text{wb}})$<br/>Coupled $[V(t), T(t)]$ | Calibrated algebraic and 2-ODE surrogates segmented by $\text{Re}$ with up to 600× speedup | JAIIO 2026 reports in `JAIIO/` |
| **`Neural`** | `Option3bNeuralODE`, `Option3cHybridODE` | SciML ODE integration | Pure MLP ODE or Physics + Transformer Regime-Attention residual | JAIIO 2026 reports in `JAIIO/` |

---

## 2. Model Family Checklist: Implemented (`[x]`) vs. Missing (`[ ]`)

### Family A: Classical & Extended Film-Theory (Spherical In-Flight Droplets)
*Papers: Langmuir (1918), Godsave (1953), Spalding (1953), Aroesty (1969), Law (1982), Abramzon & Sirignano (1989), Aggarwal & Peng (1995), Sirignano (2010), Hussain et al. (2018), Akkuş (2020), Dalla Barba et al. (2021), Xu & Xu (2025)*

- `[x]` Stationary diffusive mass loss $\dot{m} \propto r_s$ (`:d2_classic`, Langmuir 1918).
- `[x]` Spalding mass & thermal numbers $B_M, B_T$ and logarithmic Stefan flow $\ln(1+B_M)$ (`:d2_stefan`, Spalding 1953).
- `[x]` Abramzon-Sirignano blowing film factor $F(B) = (1+B)^{0.7}\frac{\ln(1+B)}{B}$, modified $\text{Sh}^*, \text{Nu}^*$, and iterative non-unitary Lewis coupling $B_T = (1+B_M)^\phi - 1$ (`solve_as1989_transfer`).
- `[x]` Frössling/Ranz-Marshall and Clift-Grace-Weber convective correlations (`nusselt_froessling`, `clift_correlation`).
- `[x]` Effective thermal conductivity $\chi_l = 1.86 + 0.86\tanh[2.245\log_{10}(\text{Pe}_l/30)]$ and blowing friction drag $C_F = \frac{12.69}{\text{Re}^{2/3}(1+B_M)}$ (`internal_circulation_enhancement`).
- `[x]` 1/3 Reference State Rule (Abramzon & Sirignano 1989 Eq. 4): `reference_state_one_third` and `reference_state = :one_third` in `DropletParams`.
- `[x]` Wilke Gas-Mixture Rules (Abramzon & Sirignano 1989 Appendix): `wilke_mixture_viscosity`, `wilke_mixture_conductivity`, `gas_mixture_properties`, and `use_wilke_mixtures = true` in `DropletParams`.
- `[x]` Revised $D^2$-Law (`:d2_revised_dallabarba`) & Ventilated Wet-Bulb Equilibrium: `ventilated_wet_bulb_temperature` and `calc_d2_revised_rate` (Dalla Barba et al. 2021, Aroesty 1969, Kinzer & Gunn 1951).
- `[x]` Non-Spherical Spheroid Deformation Factor: `spheroid_deformation_evaporation_ratio` and `spheroid_epsilon` in `DropletParams` (Tonini & Cossali 2015, Sazhin 2021).
- `[ ]` **Dynamic Drag & Deceleration ODEs (Abramzon 1989 Eqs. 1–2, Miesse 1954, Aggarwal 1995)**: Support dynamic droplet deceleration $du/dt = \frac{3}{8}\frac{\rho_g}{\rho_l}\frac{C_D}{r_s}|U_\infty - u|(U_\infty - u)$ in addition to instantaneous terminal velocity.
- `[x]` **Dimensionless Gas Conduction Number $N_{hc}$ (Xu & Xu 2025)**: Diagnostic evaluating gas-phase conduction vs. latent cooling (`xu_gas_conduction_number`).

---

### Family B: Internal Liquid Heating & Temperature Profile Models
*Papers: Sazhin (2006, 2017, 2021), Zhou et al. (2017)*

- `[x]` Uniform internal temperature / infinite conductivity (`:uniforme`, $b_q = 0$).
- `[x]` Parabolic temperature profile (`:parabolico`, $b_q = 5\chi_l$) and cubic coefficient (`:cubico`, $b_q = 11.5$).
- `[x]` 1D Moving-boundary radial conduction PDE with immobilized Landau coordinate $\xi = r/R(t)$ (`Evaporacion.RadialPDE`).
- `[x]` Zhou et al. (2017) Closed-Form 3rd-Order Profile (Zhou 2017 Eqs. 11–18): `reconstruct_zhou_cubic_profile`, `:zhou_cubic` in `DropletParams`, moving-boundary $(T_s - \bar{T})$ term in `dT/dt`, and Hill's vortex stagnation point at $r = r_s/\sqrt{2}$.
- `[x]` Angelilli et al. (2024) Liquid Thermal Dilation & Variable-$c_{p,l}$ Factor: `liquid_thermal_expansion_water`, `thermal_dilation = true`, `variable_cp_factor = true`.
- `[x]` Sazhin's Analytical Eigenvalue Series (Sazhin 2006, 2017, 2021): `sazhin_eigenvalue_temperature_profile` evaluating exact analytical series with bracketed root-finding for $\tan\lambda_n = -\lambda_n/(\text{Bi}-1)$.
- `[ ]` **Thermal Radiation Absorption (Sazhin 2006 Section 3.5)**: Droplet radiative heating $\dot{q}_{\text{rad}} = \pi R_d^2 \bar{Q}_{\text{abs}} \sigma_{\text{SB}}(T_{\text{rad}}^4 - T_s^4)$ with $\bar{Q}_{\text{abs}} = 1 - e^{-a_r R_d}$.

---

### Family C: Agricultural Spray Drift, Atmospheric Sedimentation & Ballistics
*Papers: Kinzer & Gunn (1951), Miesse (1954), Holterman (2003), Renaudo et al. (2022, 2023, 2024, 2025), Okada et al. (2021), Chiu / Sirignano (2010)*

- `[x]` Stull (2011) psychrometric wet-bulb temperature (`stull_twb`).
- `[x]` Dynamic free-fall terminal velocity in air (`calc_terminal_velocity`).
- `[x]` Reproduction of Renaudo et al. (2025) CLICAP case studies in [`clicap/reproduce_paper_cases.jl`](file:///d:/programacion/evaporacion/clicap/reproduce_paper_cases.jl).
- `[x]` Holterman (2003) $(a, b)$ Weather Parameterization & Relaxation Time $\tau_c$: `holterman_evaporation_rate` ($dD/dt = -\frac{a}{D}(1 + b\sqrt{Dv})$) and `holterman_relaxation_time` ($\tau_c$).
- `[x]` **In-Flight Deposition & Extinction Equations (Renaudo et al. 2022, 2023)**: Analytical functions for evaporation parameter $k = \frac{1.25\times 10^9}{T_c - T_{bh}}$, cut-off extinction diameter $d_{\text{min}} = [\frac{54\mu_g H}{k(\rho_d-\rho_g)g}]^{1/4}$ (`renaudo_extinction_diameter`), in-flight soil deposition time $t_{\text{dep}}(d_0, H)$ (`renaudo_deposition_time`), and landing diameter $d_{\text{dep}}(d_0) = d_0 (1 - (d_{\text{min}}/d_0)^4)^{1/6}$ (`renaudo_landing_diameter`).
- `[x]` **Chiu's Group Evaporation Number (Sirignano 2010)**: Cloud/spray collective evaporation parameter $G = 3(1 + 0.276\text{Re}^{1/2}\text{Sc}^{1/3})\text{Le} N_{\text{tot}}^{2/3}(R_d/R_c)$ and regime classifier (`chiu_group_evaporation_number`, `chiu_evaporation_regime`).
- `[ ]` **Lagrangian Atmospheric Turbulent Dispersion (Renaudo 2023, 2024)**: Random-walk eddy diffusivity trajectory tracking in turbulent crosswinds.

---

### Family D: Diffusion-Limited Sessile Droplet Evaporation Modes
*Papers: Picknett & Bexon (1977), Hu & Larson (2002), Dunn et al. (2008), Stauber et al. (2015), Saxton et al. (2016), Ong (2022), Teng et al. (2023), Sourais et al. (2023), Wilson & D'Ambrosio (2023)*

- `[x]` Exact toroidal Mehler-Fock integral $g(\theta)$ (`g_theta_exact`) and fast precomputed interpolation table (`g_theta_fast`).
- `[x]` Hu-Larson (2002) and Hu (2014) empirical approximations (`g_theta_hu_larson`, `g_theta_hu_2014`).
- `[x]` Analytical lifetimes and ODE solvers for all 4 canonical modes: `:CR`, `:CA`, `:SS`, and `:SJ` (`lifetime_cr/ca/ss/sj`, `solve_sessile`).
- `[x]` Local evaporative flux profiles: `local_flux_popov_series` and `local_flux_deegan` ($J(r) \propto (1-r/R)^{\lambda-1}$).
- `[x]` **Picknett & Bexon (1977) Capacitance Polynomials**: Two-branch polynomial for equiconvex lens capacitance $C/r(\theta)$ (`picknett_capacitance_ratio`).
- `[x]` **Thin-Droplet Closed-Form Lifetimes (Wilson & D'Ambrosio 2023 Section 3.3, Dunn 2008)**: Shallow-limit ($\theta_0 \ll 1$) formulas $t_{\text{CR,thin}} = \frac{\pi\rho\theta_0 R_0^2}{16 D \Delta c}$, $t_{\text{CA,thin}} = 1.5 t_{\text{CR,thin}}$, $t_{\text{SS,thin}}$, $t_{\text{SJ,thin}}$ (`thin_droplet_lifetime_cr/ca/ss/sj`).
- `[x]` **Multi-Droplet Vapor Shielding Factor (Wilson & D'Ambrosio 2023 Eq. 38)**: Lifetime extension factor $\eta(B/R_0) = 1 + \frac{2}{\pi}\arcsin(R_0/B)$ for droplet pairs and multi-droplet arrays (`droplet_pair_shielding_factor`, `droplet_array_shielding_factor`).
- `[ ]` **Dynamic Pinning/Depinning Force Balance on Heterogeneous Substrates (Sourais et al. 2023)**: Contact line motion governed by capillary depinning force threshold $\gamma(\cos\theta - \cos\theta_{\text{eq}})$.

---

### Family E: Sessile Internal Flow, Substrate Thermal Conduction & Marangoni Effects
*Papers: Niazmand et al. (1994), Larson (2014), O'Brien (2018), Gelderblom et al. (2022), Pan & Wang (2022), Jain et al. (2024), Rocha et al. (2024), Dekker et al. (2025)*

- `[x]` Pinned (CR) lubrication velocity fields $\bar{u}(r)$, $u(r, z)$, $w(r, z)$ and coffee-ring mass deposit scaling (`internal_deposition.jl`, Wilson 2023 Eqs. 24–26).
- `[x]` Substrate conductive cooling $\Delta c_{\text{eff}}$ and temperature drop estimate (`SessileThermalParams`, `cooling_effective_delta_c`).
- `[x]` Spherical thermal Marangoni flow factor and surfactant adsorption kinetics in `Evaporacion.RadialPDE`.
- `[x]` **Unpinned (CA) Radial Velocity Field (Gelderblom et al. 2022 Eq. 3.2.2)**: Lubrication velocity profile $\bar{u}_{\text{CCA}}(r, t)$ where radial speed vanishes at the receding contact line (`lubrication_velocity_cca`).
- `[x]` **Ristenpart Thermal Inversion Criterion & Marangoni Direction (Wilson 2023 Section 3.4, Gelderblom 2022)**: Diagnostic evaluating $k_R = k_{\text{substrate}} / k_{\text{liquid}}$, thermal Marangoni number $\text{Ma}_T$, solutal Marangoni number $\text{Ma}_S$, and the $\text{Ma}/\text{Ra}$ buoyancy ratio (`ristenpart_thermal_ratio`, `thermal_marangoni_number`, `solutal_marangoni_number`, `marangoni_rayleigh_ratio`).
- `[x]` **Asymmetric Gravitational Sagging on Inclined Substrates (Pan & Wang 2022)**: Contact angle distortion ($\theta_{\text{front}} > \theta_0 > \theta_{\text{rear}}$) as a function of Bond number $\text{Bo}$ and inclination angle $\alpha$ (`inclined_substrate_sagging`).

---

### Family F: Multicomponent Droplets & Non-Ideal Mixtures
*Papers: Yao (2005), Longest & Kleinstreuer (2005), Tonini & Cossali (2015), Ordoubadi et al. (2019), Angelilli et al. (2024)*

- `[x]` Single-component water evaporation with non-volatile solute Raoult mole-fraction activity $x_w$ in `RadialPDE`.
- `[x]` **Tonini & Cossali (2015) Analytical Stefan-Maxwell Matrix Exponential Solver**: Exact matrix solution $e^{\mathbf{A}}(\mathbf{A}\mathbf{W}_\infty + \mathbf{B}) - (\mathbf{A}\mathbf{W}_s + \mathbf{B}) = \mathbf{0}$ for coupled multicomponent diffusion without effective Fickian diffusivity approximations (`stefan_maxwell_matrices`, `stefan_maxwell_residual`, `solve_stefan_maxwell_exact`).
- `[x]` **Multicomponent Species Conservation ODEs & Liquid Activity (Yao 2005, Longest 2005, Ordoubadi 2019)**: Tracking multi-species liquid mass fractions $dm_i/dt = -\dot{m}_i$ with Ordoubadi Stefan-Fuchs radial advection fluxes $\dot{m}_t$ and individual solvent mass rates (`OrdoubadiStefanFuchsFluxes`, `solve_multicomponent_droplet`).

---

### Family G: Kinetic Non-Equilibrium & Knudsen Layer Models
*Papers: Sazhin (2006), Vlasov (2021), Zhao & Nadal (2023), Wilson & D'Ambrosio (2023), Akkuş (2020)*

- `[x]` Continuum diffusion-limited boundary condition ($Y_s = Y_{\text{sat}}(T_s)$, valid for $\text{Kn} \ll 1$).
- `[x]` **Mixed Kinetic-Diffusive Interface Jump (Vlasov 2021, Zhao & Nadal 2023, Sazhin 2006)**: Hertz-Knudsen-Schrage boundary resistance for micro/nanodroplets ($d \lesssim 10\ \mu\text{m}$) as a function of Knudsen number $\text{Kn} = \lambda_{\text{mfp}}/R$ and accommodation coefficient $\beta_{\text{evap}}$ (`vlasov_kinetic_rate`, `vlasov_kinetic_length`, `vlasov_droplet_radius`, `vlasov_lifetime_extension`, `zhao_peclet_number`, `zhao_lifetime_extension`, `solve_kinetic_droplet`).
- `[x]` **Kelvin Curvature Elevation & Stefan Blowing (Akkuş 2020)**: Elevated equilibrium vapor pressure $p_v(r) = p_{\text{sat}}\exp(\frac{2\sigma M_v}{\rho_l R_g T_s r})$ and high-mass-flux blowing enhancement $\Theta_{\text{Stefan}} = \frac{1}{1 - Y_{vs}}$ (`akkus_kelvin_vapor_ratio`, `akkus_stefan_blowing_factor`).

---

### Family H: Particle Suspensions, Nanofluids & Salt Crusts
*Papers: Gerken et al. (2014), Sadek et al. (2015), Maruyama & Hasegawa (2020), Ochowiak et al. (2022), Ordoubadi et al. (2019), Vehring et al. (2007)*

- `[x]` Ideal Raoult solute activity in `RadialPDE`.
- `[x]` **Modified $D^2$-Law with Insoluble Solid Fraction & Fractal Crusts (Gerken 2014, Ochowiak 2022)**: Fractal agglomerate packing limit $\phi_{\text{max}} = (R_g/r_{np})^{d_f-3}\frac{\pi}{6}$ (`fractal_packing_limit`), locking final dry-particle diameter $d_{\text{lock}} = d_0 \phi_{s,0}^{1/3}$ (`crust_locking_diameter`).
- `[x]` **Electrolyte Water Activity & Shell Formation Peclet Number (Maruyama 2020, Vehring 2007)**: Raoult-van 't Hoff vapor pressure depression $p_{vs} = p_{\text{sat}}(1 - i Z_s)$, time-varying evaporation constant $\beta(t)$, salt crystallization time $t_p$, dry salt core diameter $d_{\text{salt}} = d_0 (Y_{s0}\rho_{\text{sol}}/\rho_{\text{salt}})^{1/3}$ (`maruyama_vapor_pressure`, `maruyama_evaporation_rate_factor`, `solve_saline_droplet`), and Vehring solute Péclet number $\text{Pe}_j = \frac{\kappa}{8 D_j}$ with surface enrichment $E_j$ (`solute_peclet_number`, `solute_surface_enrichment`).

---

### Family I: Purely Experimental Benchmark Data
*Papers: Cioulachtjian et al. (2010), Xu et al. (2011), Devlin et al. (2015), Misyura (2017), Moon et al. (2021), Radmilovic-Radjenovic et al. (2021), Liu et al. (2023), Sun et al. (2025)*

- `[x]` Identified and cataloged in [`docs/PAPERS_CATALOG.md`](file:///d:/programacion/evaporacion/docs/PAPERS_CATALOG.md) as purely experimental reference datasets (e.g. acoustic levitation in moist air, particle deposition in sessile/pendant drops, waxy/hairy leaf adjuvant spread, heated surface interferometry).
