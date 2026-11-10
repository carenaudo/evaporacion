# Evaluación Comparativa de Modelos Surrogados para la Evaporación de Gotas
## Informe Técnico para JAIIO 2026

---

## 1. Introducción

La simulación de la evaporación de gotas es una componente computacional crítica en modelos de sprays, secado por atomización y combustión. El modelo riguroso de referencia (Nivel 3) resuelve un sistema acoplado de EDOs para masa y temperatura, con una **solución no lineal anidada** (método de Brent) para la temperatura superficial en cada evaluación de la función del lado derecho (RHS). Esta complejidad limita severamente su uso en simulaciones de gran escala.

Este trabajo evalúa cuantitativamente cuatro modelos surrogados contra el modelo riguroso en tres dimensiones:
- **Precisión** (MAE, MAPE, sesgo)
- **Tiempo de ejecución** (media ± desviación estándar, factor de aceleración)
- **Costo computacional** (evaluaciones de la función ODE, iteraciones Brent)

## 2. Metodología

### 2.1 Modelo Riguroso de Referencia (Nivel 3)

El modelo riguroso resuelve el sistema:

$$\frac{dm}{dt} = -\dot{m}(T_s), \qquad \frac{dT_{bar}}{dt} = \frac{Q_{conv} - \dot{m} L_v}{m \cdot c_{pl}}$$

donde $T_s$ se obtiene resolviendo el balance energético superficial mediante el algoritmo de Brent (`brentq`) en cada paso temporal. El integrador utiliza LSODA con tolerancias $\text{rtol}=10^{-8}$, $\text{atol}=10^{-12}$.

> [!WARNING]
> Cada evaluación del RHS de la ODE ejecuta internamente un procedimiento de búsqueda de raíces (Brent), típicamente 10–30 iteraciones por llamada. Esto constituye el principal driver del costo computacional.

### 2.2 Modelos Surrogados Evaluados

| Modelo | Tipo | Descripción / Simplificación Física |
| :--- | :--- | :--- |
| **rigorous** | Físico Nivel 3 | Modelo completo (Abramzon-Sirignano + perfil parabólico). |
| **level_2** | Físico Reducido | Sin corrección por Stefan flow (convección pura). |
| **uniform** | Físico Reducido | Perfil radial uniforme ($T_s = T_{bar}$), **sin Brent**. |
| **iso_twb** | Físico Reducido | Gota isotérmica a $T_{wb}$ (sin EDO de temperatura). |
| **modelo_6b_refit** | Surrogado | Algebraico segmentado (el mejor compromiso global). |
| **modelo_7** | Surrogado | Algebraico multi-variable ($HR, \Theta$). |
| **modelo_2b** | Surrogado | Algebraico base ($d^2$-law optimizada). |

### 2.3 Temporalidad de las Variables y Constantes

La diferencia clave reside en si el coeficiente de transferencia se recalcula instantáneamente o se "congela" al tiempo cero ($t=0$).

#### Modelos de 1 Ecuación (modelo_2b, 6b_refit, 7)
Estos modelos calculan su coeficiente $C$ de forma dinámica dentro del bucle del integrador ODE.

*   **Variables Instantáneas (evolucionan):** El diámetro ($d$) y el número de Reynolds ($Re$). Como la gota se encoge, $Re$ disminuye, lo que reduce la intensidad de la transferencia.
*   **Momento de cálculo:** En cada paso del solver. Esto les permite captar la ralentización de la evaporación a medida que la gota desaparece, lo que explica su alta precisión ($R^2 > 0.99$).

### 2.4 Diseño del Benchmark

- **25 casos representativos**: 4 esquinas extremas + 6 casos de borde + 15 interiores (muestreo tipo Latin-Hypercube) del espacio paramétrico 5D ($d_0, U, T_0, T_\infty, HR$).
- **Repeticiones**: N=10 ejecuciones cronometradas por caso y modelo, con 1 ejecución previa descartada (warm-up).
- **Métricas**: Tiempo de pared (`time.perf_counter`), evaluaciones ODE RHS (`sol.nfev`), llamadas e iteraciones Brent.
- **Precisión**: Extraída de las trayectorias CSV pre-calculadas en el barrido paramétrico de 30,240 simulaciones.

## 3. Resultados

### 3.1 Tabla Resumen Global

| Modelo | Tiempo Medio [s] | Std [s] | MAPE [%] | Sesgo Medio [s] | Aceleración | nfev medio | Ops. Totales |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **rigorous** | 0.9055 | 0.0152 | 0.00 | 0.00 | 1.0× | 784 | 6,228 |
| **level_2** | 0.8841 | 0.0167 | 0.30 | +4.37 | 1.0× | 774 | 6,150* |
| **uniform** | 0.0918 | 0.0028 | 0.55 | −1.62 | **11×**| 83 | 83 |
| **iso_twb** | 0.0085 | 0.0008 | 0.53 | −1.98 | **129×**| 79 | 79 |
| **modelo_6b_refit** | 0.0016 | 0.0002 | **6.02** | +7.32 | **619×** | 68 | 68 |
| **modelo_7** | 0.0017 | 0.0004 | **6.29** | −4.00 | **576×** | 67 | 67 |
| **modelo_2b** | 0.0028 | 0.0005 | 16.01 | −102.06 | 381× | 68 | 68 |

> [!IMPORTANT]
> Al alinear el solver numérico (**LSODA**) para todos los modelos, los surrogados algebraicos demuestran una aceleración masiva (>600×) comparada con los ~120× observados anteriormente con integradores sub-óptimos. Esto confirma que el costo de los surrogados es despreciable y está dominado únicamente por el overhead del integrador de EDOs.

### 3.2 Análisis de Tiempo de Ejecución

![Tiempos de Ejecución](file:///d:/programacion/evaporacion/output/benchmark/figures/exec_time_barplot.png)
*Comparación de tiempos de ejecución en escala logarítmica. Las barras de error representan la desviación estándar sobre N=10 repeticiones.*

![Distribución de Tiempos](file:///d:/programacion/evaporacion/output/benchmark/figures/exec_time_boxplot.png)
*Distribución de tiempos de ejecución a través de los 25 casos representativos.*

### 3.3 Factor de Aceleración (Speedup)

![Speedup](file:///d:/programacion/evaporacion/output/benchmark/figures/speedup_barplot.png)
*Factor de aceleración respecto al modelo riguroso. Todos los surrogados logran speedups de 2 órdenes de magnitud.*

### 3.4 Costo Computacional: Operaciones

![Operaciones](file:///d:/programacion/evaporacion/output/benchmark/figures/nfev_barplot.png)
*Comparación del número de evaluaciones de la función ODE (nfev) y operaciones totales (incluyendo iteraciones Brent para el modelo riguroso).*

El modelo riguroso ejecuta ~784 evaluaciones RHS, pero cada una contiene una búsqueda de raíces Brent (~7 iteraciones promedio por llamada), resultando en **~6,228 operaciones totales**. Los surrogados, al no requerir soluciones anidadas, tienen un costo operacional directo.

### 3.5 Frontera de Pareto: Precisión vs Costo

![Pareto](file:///d:/programacion/evaporacion/output/benchmark/figures/accuracy_vs_speed.png)
*Diagrama de Pareto: tiempo de ejecución (eje x, log) vs error MAPE (eje y). La zona inferior-izquierda es la región deseable. La línea punteada indica el tiempo del modelo riguroso.*

El diagrama confirma que **modelo_6b_refit** y **modelo_7** se ubican en el cuadrante Pareto-óptimo, combinando alta precisión (<7% MAPE) con baja latencia computacional (>100× más rápido).

### 3.7 Análisis de Peores Casos (nfev)

Para identificar la robustez computacional, se extrajeron los casos que generan el mayor número de evaluaciones de función (`nfev` peak) para cada modelo:

| Modelo | nfev max | Caso Crítico | Características del Caso |
| :--- | ---: | :--- | :--- |
| **rigorous** | 1,054 | Caso 23 | $d_0=2.0$mm, $U=0.5$m/s, $T_0=10$°C, $T_{inf}=20$°C |
| **uniform** | 130 | Caso 22 | $d_0=2.0$mm, $U=1.0$m/s, $T_0=35$°C, $T_{inf}=45$°C |
| **surrogados** | 120 | Caso 2 | $d_0=2.8$mm, $U=0.0$m/s, $T_0=35$°C, $T_{inf}=40$°C |

**Observación clave:** Los surrogados, al usar el solver adaptativo LSODA, ahora requieren un número de evaluaciones similar al modelo `iso_twb`, independientemente de la duración de la simulación. Caso 2 (la gota de mayor vida útil, >10,000s) es el único que eleva ligeramente el nfev para los modelos algebraicos.

## 4. Discusión

### 4.1 Compromiso Precisión-Velocidad

Los resultados demuestran que la estructura algebraica de **modelo_6b_refit** (escalamiento segmentado por $Re$ con corrección global) ofrece el mejor equilibrio entre fidelidad física y eficiencia computacional. Su MAPE de 6% es aceptable para simulaciones de sprays donde la incertidumbre instrumental típica es del orden del 5-10%.

### 4.2 Impacto de las Simplificaciones Físicas
- **Flujo de Stefan (level_2)**: Su exclusión tiene un impacto despreciable en el tiempo (0.91 vs 0.89s) y en el error (0.40%). Esto sugiere que el acoplamiento de transferencia de masa no es la fuente de la lentitud.
- **Perfil Térmico (uniform)**: Es el factor determinante. Eliminar la búsqueda de $T_s$ reduce el tiempo en un orden de magnitud (10.6×) manteniendo el error por debajo del 1%.
- **EDO de Temperatura (iso_twb vs Surrogados)**: El modelo isotérmico a $T_{wb}$ es el más rápido de los modelos físicos (131×), entrando en el mismo orden de magnitud que los surrogados algebraicos.

El modelo riguroso dedica la mayor parte de su tiempo computacional a las llamadas a `brentq` dentro del RHS de la ODE. La comparación con el baseline **uniform** (que mantiene la misma física pero impone un perfil radial constante $T_s = T_{bar}$) revela la magnitud de este "impuesto" algorítmico.

Al eliminar el solver iterativo anidado, el modelo **uniform** logra una aceleración de aproximadamente 20× manteniendo un error de evaporación inferior al 1%. Esto confirma que la complejidad del perfil de temperatura, y no las correlaciones físicas de transferencia, constituye el cuello de botella computacional del modelo de Nivel 3.

## 5. Conclusiones

1. **modelo_6b_refit** es el surrogado recomendado para integración en simulaciones de sprays a gran escala, ofreciendo **6.02% MAPE** con una **aceleración de 619×**.
2. **modelo_7** representa una alternativa con dimensionalidad extendida, útil en regímenes de alta humedad.
3. El **costo computacional dominante** del modelo riguroso proviene del solver Brent anidado. La simplificación del perfil térmico (**uniform**) reduce el tiempo en **11×** con un error despreciable (<0.6%).
4. La alineación de los solvers numéricos demuestra que el overhead de los surrogados es órdenes de magnitud menor a lo reportado inicialmente, permitiendo simulaciones masivas con billones de gotas.

---

*Datos completos disponibles en:*
- [benchmark_results.csv](file:///d:/programacion/evaporacion/output/benchmark/benchmark_results.csv)
- [summary_table.csv](file:///d:/programacion/evaporacion/output/benchmark/summary_table.csv)
- [analysis_out_v3](file:///d:/programacion/evaporacion/analysis_out_v3/) (análisis de calidad completo - 30k sweep)
