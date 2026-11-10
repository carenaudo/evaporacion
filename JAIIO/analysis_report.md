# Neural ODE Model Analysis Report

## 1. Multi-Case Validation Summary
The model was validated against the physical model across 19 distinct initial conditions (d0: 0.1–2.8 mm, U: 0–3 m/s, HR: 30–70%).

| Metric | Mean Error | Max Error | Note |
| :--- | :--- | :--- | :--- |
| **Diameter (MAE)** | 0.0542 mm | 0.1520 mm | High precision (~2.5% relative error for large drops) |
| **Temperature (MAE)** | 1.0069 °C | 7.2937 °C | Good tracking; highest error in impulsive small-drop cases |
| **Evaporation Time** | 0% | 0% | Integration is stable over the full physical window |

> [!NOTE]
> The model captures the $d^2$ law dynamics and the thermal cooling/reheating cycle with remarkably few parameters (~36k).

## 2. Overfitting Analysis
Based on the training history (80 epochs):
* **Train Loss Final**: ~0.1374
* **Val Loss Final**: ~0.1387
* **Gap**: < 1%

**Conclusion**: **NO significant overfitting** is present. The validation loss tracks the training loss very closely. The current model is actually "under-capacity" or "under-trained" rather than overfitted. We can afford to increase network complexity or train for more epochs (e.g., 500+) without immediate risk of divergence.

## 3. Proposed Architecture Changes
While the current MLP-based Neural ODE is functional, the following changes are recommended to reach "production-grade" precision:

### A. Physical Residual (Hybrid Model)
Instead of learning the full $dm/dt$ and $dT/dt$, the network should learn the **residual** relative to a simple $d^2$ law base.
* **Logic**: $f_{total} = f_{analytic}(d, T) + \text{NeuralNet}(d, T, \text{cond})$
* **Benefit**: Guarantees basic physical behavior even in regions with sparse data.

### B. Time-Aware Attention
The current model uses a static condition vector. Adding an **Attention Mechanism** on the condition vector relative to the current state $t/t_{end}$ would help the model "focus" on different physical regimes (e.g., the initial heat-up vs. final dry-out).

### C. Singularity-Aware Activation
The evaporation rate often scales as $1/d$. The current MLP has to "work hard" to approximate this with standard Sigmoid/SiLU activations.
* **Proposal**: Use a custom activation layer for the mass derivative that includes a $1/d$ term:
  $$\frac{dm}{dt} = -\text{Softplus}(\text{MLP}) \cdot \left( \frac{1}{d + \epsilon} \right)$$

### D. Multi-Head Condition Encoder
Instead of a single bottleneck (32-dim), use separate encoders for **Geometry** (d0, t_end) and **Environment** (Tinf, HR, U). This prevents the "Geometry" signal from being washed out by "Environment" noise.

---
**Validation Table CSV**: [validation_table.csv](file:///D:/programacion/evaporacion/saved_models/validation_table.csv)
**Full Trajectory Plots**: [validation_trajectories.png](file:///D:/programacion/evaporacion/saved_models/validation_trajectories.png)
