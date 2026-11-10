# Physics-Informed Neural Networks (PINNs) & Operator Surrogates Report (Julia)

*Generated from [train_jaiio_pinns.jl](file:///d:/programacion/evaporacion/JAIIO/train_jaiio_pinns.jl) and [evaluate_jaiio_pinns.jl](file:///d:/programacion/evaporacion/JAIIO/evaluate_jaiio_pinns.jl) (`Evaporacion.Neural.PINNSurrogates`) against the Rigorous Level 3 Julia baseline (`:brent`, Zhou et al., 2017 moving-boundary).*

---

## 1. Diagnosis of Legacy Neural Approaches & Design Principles

An audit of the legacy neural implementations (`option3a_neural_rhs.py`, `option3b_neural_ode.py`, `option3c_hybrid_attn.jl`) identified two structural flaws:
1. **Target Leakage in Legacy `Option 3b`**: The network conditioned on `log_t_end` (the ground-truth `t_evap`), preventing autonomous evaporation-time prediction, while integrating $\ln(m) \to -\infty$ introduced dry-out stiffness.
2. **Slow Runtime in `Option 3c`**: Evaluating the full iterative surface-temperature root-finder inside every neural ODE step made inference slower than the rigorous baseline.

The new Julia neural module ([pinn_surrogates.jl](file:///d:/programacion/evaporacion/julia/src/neural/pinn_surrogates.jl)) enforces **strict zero-leakage autonomy** (inputs are strictly $(d_0, U_\infty, T_0, T_\infty, \text{HR})$) across four architectures:

| Architecture | Internal ID | Formulation & Physics Constraints | Params | ODE Solver |
| :--- | :--- | :--- | ---: | :--- |
| **`PINN-1Eq`** | `pinn_c` | Smooth $C^\infty$ `modelo_7` prior $\times \exp(0.15\tanh(\mathcal{N}_\theta))$ with dimensionless $\text{Re}$ gates; trained via 8-pt Gauss-Legendre $t_{\text{evap}}$ integral + $\partial C/\partial U \ge 0$, $\partial C/\partial \text{HR} \le 0$ monotonicity penalties | ~865 | `Tsit5()` |
| **`PINN-2Eq`** | `pinn_2eq` | Coupled $[V(t), \bar{T}(t)]$ conservation with 0-iter Neural Biot $\hat{T}_s$ closure, film Sherwood/Spalding blowing, and regularized thermal inertia; trained on 24-pt collocation trajectories + energy balance & lifetime integral | ~1,507 | `AutoTsit5(Rosenbrock23)` |
| **`PI-DeepONet`** | `pi_deeponet` | Solver-free Branch-Trunk operator network (`nfev = 1`) with hard BCs ($\hat{d}(0)=d_0, \hat{d}(1)=0, \hat{\bar{T}}(0)=T_0$) and interior ODE collocation loss | ~6,350 | None (`nfev = 1`) |
| **`ResMLP-Puro`** | `resmlp` | Black-box Residual MLP predicting $\log(t_{\text{evap}}), T_{s,0}, \log(\dot{m}_0)$ without physics priors or PDE residuals (ablation control) | ~11,500 | None (`nfev = 1`) |

---

## 2. 5-Fold Cross-Validation on V2 (2,520 Rigorous Cases)

All models were cross-validated using strict case-level 5-fold splits (2,016 train / 504 validation cases per fold):

| Model | CV $t_{\text{evap}}$ MAPE [%] | CV Max ARE [%] | CV $T_{s,0}$ MAE [°C] | CV $\dot{m}_0$ MAPE [%] | CV $d(t)$ MAE [mm] | CV $\bar{T}(t)$ MAE [°C] |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| **PINN-1Eq** (`pinn_c`) | **0.34%** | 3.02% | 10.31 | 304.31% | **0.0036** | — |
| **PINN-2Eq** (`pinn_2eq`) | **0.64%** | **1.53%** | **0.29** | **5.22%** | 0.0076 | **0.024** |
| **PI-DeepONet** (`pi_deeponet`) | 0.94% | 4.48% | 0.36 | 7.90% | 0.0252 | 2.420 |
| **ResMLP-Puro** (`resmlp`) | 8.74% | 47.66% | 1.21 | 548.48% | — | — |

![PINN Loss Convergence](figures/pinn/pinn_loss_convergence.png)

---

## 3. Global Accuracy & Speedup Benchmark (25-Case CWU + Full V2 Suite)

Combined benchmark across the 25 representative cases ($N=10$ reps) and the full 2,520-case V2 dataset:

| Model | Type | Exec Time [s] | Speedup | nfev med. | CWU | V2 $t_{\text{evap}}$ MAPE [%] | V2 Max ARE [%] | V2 $\dot{m}_0$ MAPE [%] | V2 $T_{s,0}$ MAE [°C] |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **Riguroso** | Physical | 0.03798 | 1.0× | 7468 | 76539 | 0.00% | 0.00% | 0.00% | 0.00 |
| **No-Stefan** (`level_2`) | Physical Red. | 0.01821 | 2.0× | 7420 | 12532 | 0.26% | 0.80% | 1.75% | 0.02 |
| **PINN-1Eq** (`pinn_c`) | PINN (1-Eq) | 0.00056 | **69.1×** | 303 | 5890 | **0.31%** | 2.84% | 303.78% | 10.31 |
| **PINN-2Eq** (`pinn_2eq`) | PINN (2-Eq) | 0.00654 | **5.9×** | 1792 | 60449 | **0.59%** | **1.41%** | **7.85%** | **0.36** |
| **Aproximado-3** (`modelo_7`) | Algebraic | 0.00023 | 169.1× | 309 | 429 | 0.91% | 4.09% | 303.89% | 10.31 |
| **PI-DeepONet** (`pi_deeponet`) | Operator | **0.000015** | **2535.9×** | **1** | **167** | **0.92%** | 4.21% | **8.05%** | **0.36** |
| **Aproximado-2** (`modelo_6b_refit`) | Algebraic | 0.00022 | 176.8× | 309 | 307 | 1.35% | 4.11% | 306.17% | 10.31 |
| **ResMLP-Puro** (`resmlp`) | Black-Box NN | 0.000009 | 4255.9× | 1 | 207 | 9.59% | 33.49% | 25425.80% | 1.06 |

![Pareto Front](figures/pinn/pinn_pareto_front.png)

![Parity Comparison](figures/pinn/pinn_parity_comparison.png)

---

## 4. Coupled 2-Equation Trajectory Fidelity (100 Stratified V2 Cases)

Comparison of full time-series predictions $[d(t), \bar{T}(t), T_s(t)]$ against Rigorous `:brent` trajectories:

| 2-Eq Model | $d(t)$ MAE [mm] | $\bar{T}(t)$ MAE [°C] | $T_s(t)$ MAE [°C] |
| :--- | ---: | ---: | ---: |
| **PINN-2Eq** (`pinn_2eq`) | **0.0068** | **0.024** | **0.033** |
| **PI-DeepONet** (`pi_deeponet`) | 0.0236 | 2.338 | 2.155 |

![2-Eq Trajectories](figures/pinn/pinn_trajectories_2eq.png)

---

## 5. Out-of-Distribution Generalization on V1 (640 Unseen Cases)

All models trained exclusively on V2 ($T_\infty \in [15, 50]\,^\circ\text{C}$, $\text{HR} \in [0.3, 0.7]$) were evaluated out-of-distribution on V1 ($T_\infty$ up to **$65\,^\circ\text{C}$**, $\text{HR} \in [0.1, 0.9]$):

| Model | V2 In-Dist. MAPE [%] | V1 OOD MAPE [%] | V1 OOD Max ARE [%] | V1 OOD $\dot{m}_0$ MAPE [%] | V1 OOD $T_{s,0}$ MAE [°C] | V1 OOD $R^2$ |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| **No-Stefan** (`level_2`) | 0.26% | 0.65% | 2.11% | 3.48% | 0.11 | 0.99995 |
| **PINN-2Eq** (`pinn_2eq`) | **0.59%** | **0.89%** | **2.38%** | **6.96%** | **0.96** | **0.99996** |
| **PINN-1Eq** (`pinn_c`) | 0.31% | 4.09% | 27.03% | 286.72% | 17.49 | 0.98863 |
| **Aproximado-2** (`modelo_6b_refit`) | 1.35% | 5.19% | 13.16% | 298.89% | 17.49 | 0.99589 |
| **Aproximado-3** (`modelo_7`) | 0.91% | 5.91% | 34.49% | 294.49% | 17.49 | 0.98040 |
| **PI-DeepONet** (`pi_deeponet`) | 0.92% | 6.08% | 68.24% | 7.15% | 0.96 | 0.93843 |
| **ResMLP-Puro** (`resmlp`) | 9.59% | 18.33% | 77.88% | 133.12% | 2.01 | 0.85888 |

![V1 Out-of-Distribution Generalization](figures/pinn/pinn_v1_extrapolation.png)

---

## 6. Key Engineering Conclusions

1. **Best Overall Coupled Physics Surrogate (`PINN-2Eq`)**: Achieves **0.59% V2 MAPE** (`1.41%` max error), **0.024 °C** trajectory temperature accuracy, **7.85%** $\dot{m}_0$ MAPE (reducing initial flux error by **~40×** compared to `modelo_7`), and **0.89% V1 OOD MAPE** (`2.38%` max error) when extrapolating to $T_\infty = 65^\circ\text{C}$ and $\text{HR}=10\text{–}90\%$.
2. **Highest-Accuracy 1-Equation ODE Surrogate (`PINN-1Eq`)**: Cuts `modelo_7`'s V2 MAPE by **3×** (from `0.91%` to **`0.31%`**) at **`69×` speedup** over the rigorous solver.
3. **Fastest Operator Surrogate (`PI-DeepONet`)**: Matches `modelo_7` accuracy on V2 (**0.92% MAPE**) with zero ODE solver steps (`nfev = 1`, **15 µs** per simulation, **2,536× speedup**, **167 CWU**), while improving initial flux accuracy from `304%` to **`8.05%`**.
4. **Value of Physics Constraints (`PINN` vs `ResMLP-Puro`)**: Pure black-box `ResMLP-Puro` degrades to `9.59%` in-distribution and `18.33%` out-of-distribution (`25,426%` initial flux error), confirming that embedding thermodynamic conservation laws and Biot/Spalding structure is essential for surrogate fidelity and generalization.

---

- **Core Implementation**: [pinn_surrogates.jl](file:///d:/programacion/evaporacion/julia/src/neural/pinn_surrogates.jl)
- **Training & CV Script**: [train_jaiio_pinns.jl](file:///d:/programacion/evaporacion/JAIIO/train_jaiio_pinns.jl)
- **Benchmark & Evaluation Script**: [evaluate_jaiio_pinns.jl](file:///d:/programacion/evaporacion/JAIIO/evaluate_jaiio_pinns.jl)
- **CSV Metrics**: [pinn_benchmark_summary.csv](file:///d:/programacion/evaporacion/JAIIO/results/pinn_benchmark_summary.csv), [pinn_cv_metrics.csv](file:///d:/programacion/evaporacion/JAIIO/results/pinn_cv_metrics.csv), [pinn_v2_all_metrics.csv](file:///d:/programacion/evaporacion/JAIIO/results/pinn_v2_all_metrics.csv), [pinn_v1_all_metrics.csv](file:///d:/programacion/evaporacion/JAIIO/results/pinn_v1_all_metrics.csv)
