# Modelo Riguroso de Referencia (Nivel 3)

Este documento detalla las ecuaciones y el procedimiento computacional del **Modelo Riguroso de Referencia**, basado en la formulación de **Abramzon & Sirignano (1989)**. Este modelo se utiliza como estándar de "ground truth" para evaluar la precisión de los modelos algebraicos surrogados.

## 1. Sistema de Ecuaciones Diferenciales Ordinarias (EDO)

El estado de la gota se describe mediante su masa ($m$) y su temperatura media ponderada ($\bar{T}_l$). El sistema evoluciona según:

$$ \frac{dm}{dt} = -\dot{m} $$
$$ \frac{d\bar{T}_l}{dt} = \frac{Q_{conv} - \dot{m} L_v}{m \cdot c_{pl}} $$

Donde:
- $\dot{m}$: Tasa de evaporación másica [kg/s].
- $Q_{conv}$: Flujo de calor convectivo desde el gas hacia la superficie [W].
- $L_v$: Entalpía de vaporización evalua a la temperatura superficial $T_s$ [J/kg].
- $c_{pl}$: Calor específico del líquido [J/kg·K].

### 2.1 Números Adimensionales y Transferencia Base
Antes de calcular los números de transferencia, se definen los parámetros adimensionales del flujo gaseoso:

1. **Número de Reynolds ($Re$)**: Relaciona fuerzas inerciales y viscosas.
   $$ Re = \frac{\rho_g \cdot |U| \cdot d}{\mu_g} $$
   Donde $|U|$ es la velocidad relativa entre la gota y el gas.
2. **Número de Prandtl ($Pr$)**: Relaciona la difusividad de momento y la térmica.
   $$ Pr = \frac{\mu_g \cdot c_{pg}}{k_g} $$
3. **Número de Schmidt ($Sc$)**: Relaciona la difusividad de momento y la de masa.
   $$ Sc = \frac{\mu_g}{\rho_g \cdot D_{ab}} $$

Con estos parámetros, se calculan los números de Nusselt ($Nu_0$) y Sherwood ($Sh_0$) base (sin corregir por soplido) mediante la correlación de Ranz-Marshall:
$$ Nu_0 = 2 + 0.552 Re^{1/2} Pr^{1/3} $$
$$ Sh_0 = 2 + 0.552 Re^{1/2} Sc^{1/3} $$

### 2.2 Factor de Soplido y Números de Spalding
Los números de Spalding para masa ($B_M$) y calor ($B_T$) cuantifican la fuerza impulsora para la transferencia:

$$ B_M = \frac{Y_s - Y_\infty}{1 - Y_s} $$

Donde:
- $Y_s$: Fracción másica de vapor en la superficie de la gota. Se asume saturación local:
  $$ Y_s = \frac{X_s \cdot M_v}{X_s \cdot M_v + (1 - X_s) \cdot M_g} $$
  siendo $X_s = P_{sat}(T_s) / P_{atm}$ la fracción molar en superficie.
- $Y_\infty$: Fracción másica de vapor en el ambiente infinito:
  $$ Y_\infty = \frac{X_\infty \cdot M_v}{X_\infty \cdot M_v + (1 - X_\infty) \cdot M_g} $$
  siendo $X_\infty = (HR \cdot P_{sat}(T_\infty)) / P_{atm}$ la fracción molar en el infinito.
- $M_v, M_g$: Masas molares del vapor (agua) y del gas portador (aire).
- $HR$: Humedad relativa del ambiente [0-1].

El número de Spalding para calor ($B_T$) se obtiene del acoplamiento:
$$ B_T = \left( 1 + B_M \right)^\phi - 1 $$

El factor de soplido universal $F(B)$ se define como:
$$ F(B) = (1 + B)^{0.7} \frac{\ln(1+B)}{B} $$

Los números de transferencia corregidos son:
$$ Sh^* = 2 + \frac{Sh_0 - 2}{F(B_M)} $$
$$ Nu^* = 2 + \frac{Nu_0 - 2}{F(B_T)} $$

### 2.3 Tasa de Evaporación y Transferencia de Calor
La tasa de evaporación másica instantánea ($\dot{m}$) y el flujo de calor convectivo ($Q_{conv}$) se calculan como:

$$ \dot{m} = \pi d \rho_g D_{ab} Sh^* \ln(1 + B_M) $$
$$ Q_{conv} = \pi d k_g Nu^* (T_\infty - T_s) \frac{\ln(1 + B_T)}{B_T} $$

Donde:
- $d$: Diámetro de la gota [m].
- $\rho_g, k_g$: Densidad y conductividad del gas en la película [kg/m³, W/m·K].
- $D_{ab}$: Difusividad binaria en la película [m²/s].
- $B_M, B_T$: Números de Spalding para masa y calor.

### 2.4 Acoplamiento Calor-Masa
La relación $\phi$ acopla ambos procesos:
$$ \phi = \frac{c_{pv}}{c_{pg}} \frac{Sh^*}{Nu^*} \frac{1}{Le} $$
Donde $Le = Sc/Pr$ es el número de Lewis. Debido a que $Nu^*$ depende de $B_T$ y $B_T$ depende de $\phi$ (que a su vez depende de $Nu^*$), este sistema se resuelve iterativamente para cada evaluación del RHS.

## 3. Resolución de la Temperatura Superficial ($T_s$)

La temperatura superficial $T_s$ no es una variable de estado directa, sino que se obtiene del balance energético en la interfase gas-líquido:

$$ Q_{conv} = \dot{m} L_v + Q_l $$

Donde $Q_l$ es el flujo de calor que penetra hacia el núcleo de la gota.

### 3.1 Perfil de Temperatura Interno (Modelo Parabólico)
El modelo asume que el gradiente térmico dentro de la gota sigue un perfil radial parabólico:
$$ T(r, t) = T_c(t) + [T_s(t) - T_c(t)] \left( \frac{r}{R} \right)^2 $$
Donde $T_c$ es la temperatura en el centro. Esta suposición permite relacionar la temperatura media ($\bar{T}_l$) con la temperatura superficial ($T_s$):
$$ \bar{T}_l = T_c + \frac{3}{5}(T_s - T_c) \implies (T_s - T_c) = \frac{5}{3}(T_s - \bar{T}_l) $$

El gradiente en la superficie ($\partial T / \partial r|_{r=R}$) resulta en un flujo de conducción:
$$ Q_{cond} = 4\pi R^2 k_l \left( \frac{\partial T}{\partial r} \right)_{r=R} = \frac{40}{3} \pi R k_l (T_s - \bar{T}_l) $$

### 3.2 Mejora por Circulación Interna ($\chi_l$)
Para gotas en movimiento, la fricción gaseosa induce una circulación interna (vórtice de Hill) que acelera la transferencia de calor. Se define una conductividad térmica efectiva $k_{eff} = \chi_l k_l$, donde $\chi_l$ se calcula iterativamente en función del número de Peclet líquido ($Pe_l$):
$$ Pe_l = \frac{U_s \cdot d}{\alpha_l} $$
Donde $\alpha_l = k_l / (\rho_l \cdot c_{pl})$ es la difusividad térmica del líquido y $U_s$ es la velocidad de la superficie del vórtice inducida por el gas:
$$ U_s = \frac{\rho_g \cdot C_F \cdot |U| \cdot d}{32 \cdot \mu_l} $$
El coeficiente de fricción corregido ($C_F$) se calcula según la correlación para gotas con soplido (blowing):
$$ C_F = \frac{12.69 \cdot Re^{-1/3}}{1 + B_M} $$
(Válido para $Re > 10^{-3}$, de lo contrario $C_F = 0$). El factor de incremento $\chi_l$ es:
$$ \chi_l = 1.86 + 0.86 \tanh \left[ 2.245 \log_{10} \left( \frac{Pe_l}{30} \right) \right] $$

El flujo de calor total al líquido se expresa entonces como:
$$ Q_l = b_q \cdot \pi R (\chi_l k_l) (T_s - \bar{T}_l) $$
Donde $b_q = 20/3$ (o un valor similar dependiendo de la aproximación del perfil, típicamente simplificado a 5.0 o 6.67 en implementaciones numéricas). En este código se utiliza un coeficiente de transferencia interna $b_q$ derivado del perfil.

### 3.3 Algoritmo de Brent
En cada evaluación del RHS de la EDO, se busca la raíz $T_s$ del residuo:
$$ \text{Residual}(T_s) = Q_{conv}(T_s) - \dot{m}(T_s) L_v(T_s) - Q_l(T_s, \bar{T}_l) = 0 $$

### 3.4 Perfiles Alternativos (Cúbico)
Además del perfil parabólico, el modelo permite seleccionar un **perfil cúbico** para la distribución interna de temperaturas. En este caso, el coeficiente de transferencia interna $b_q$ se ajusta a **11.5**, lo que modela un gradiente térmico más pronunciado en la superficie, útil para capturar la rápida penetración de calor durante los instantes iniciales de la simulación.

## 4. Propiedades Físicas y Temperaturas de Evaluación

Para asegurar la precisión del modelo, las propiedades deben evaluarse en la temperatura local correspondiente al fenómeno físico (gas, superficie o interior líquido).

### 4.1 Resumen de Evaluación por Fase

| Propiedad | Símbolo | Fase | Temperatura de Evaluación |
| :--- | :---: | :--- | :--- |
| **Densidad Gas** | $\rho_g$ | Gas | Temperatura de Película ($T_{film}$) |
| **Viscosidad Dinámica** | $\mu_g$ | Gas | Temperatura de Película ($T_{film}$) |
| **Conductividad Térmica** | $k_g$ | Gas | Temperatura de Película ($T_{film}$) |
| **Calor Específico Gas** | $c_{pg}$ | Gas | Temperatura de Película ($T_{film}$) |
| **Calor Específico Vapor** | $c_{pv}$ | Gas | Temperatura de Película ($T_{film}$) |
| **Difusividad Binaria** | $D_{ab}$ | Gas | Temperatura de Película ($T_{film}$) |
| **Presión de Vapor** | $P_{sat}$ | Interfase | Temperatura Superficial ($T_s$) |
| **Entalpía Vaporización** | $L_v$ | Interfase | Temperatura Superficial ($T_s$) |
| **Densidad Líquido** | $\rho_l$ | Líquido | Temperatura Media ($\bar{T}_l$) |
| **Calor Especifico Líquido** | $c_{pl}$ | Líquido | Temperatura Media ($\bar{T}_l$) |
| **Conductividad Líquido** | $k_l$ | Líquido | Temperatura Media ($\bar{T}_l$) |
| **Viscosidad Líquido** | $\mu_l$ | Líquido | Temperatura Media ($\bar{T}_l$) |
| **Presión Atmosférica** | $P_{atm}$ | - | Constante (101325 Pa) |
| **Masas Molares** | $M_v, M_g$ | - | Constante (Agua, Aire) |

### 4.2 Definiciones de Temperatura

*   **Temperatura de Película ($T_{film}$)**: Es el promedio aritmético entre la superficie de la gota y el gas infinito. Se utiliza para captar las propiedades de transporte en la capa límite.
    $$ T_{film} = \frac{T_s + T_\infty}{2} $$
*   **Temperatura Media del Líquido ($\bar{T}_l$)**: Representa la temperatura uniforme equivalente de la masa de la gota. Se utiliza para las propiedades térmicas y mecánicas del líquido.
*   **Temperatura de Referencia para Vaporización**: Los efectos terminales (fracción másica en superficie $Y_s$ y calor latente $L_v$) se evalúan estrictamente a la temperatura instantánea de la interfase, $T_s$.

## 5. Complejidad Computacional

Este modelo es computacionalmente costoso debido a:
1.  **Iteración de Brent**: Búsqueda de raíces para $T_s$ en cada llamada al solver.
2.  **Iteración de AS89**: Resolución del acoplamiento $B_T - Nu^*$ dentro de cada paso de Brent.
3.  **Propiedades variables**: Recálculo constante de densidades, conductividades y difusividades.

Esto resulta en aproximadamente 10–30 iteraciones internas por cada paso del integrador ODE, lo que motiva la búsqueda de modelos surrogados eficientes.
