# Abramzon & Sirignano (1989) – Complete Extracted Equations
*(Including auxiliary Wilke/mixture equations; excluding fuel-specific polynomial fits)*

## Nomenclature (from Abramzon & Sirignano, 1989)

| Symbol | Meaning |
|---|---|
| `A_f` | Averaging parameter |
| `B_M` | Mass transfer number |
| `B_T` | Heat transfer number |
| `C_D` | Droplet drag coefficient |
| `C_f` | Friction drag coefficient |
| `c_pF` | Fuel vapor specific heat |
| `D` | Vapor/air binary diffusion coefficient |
| `F` | Relative change of film thickness |
| `i` | Iteration number |
| `k` | Thermal conductivity |
| `Le` | Lewis number, `k/(ρ D c_p)` |
| `ṁ` | Vaporization rate `[g s^-1]` |
| `m_0` | Initial droplet mass `[g]` |
| `Nu` | Nusselt number |
| `P` | Pressure |
| `Pe` | Peclet number, `Re Pr` |
| `Pr` | Prandtl number |
| `Q_L` | Heat transferred into the droplet `[cal s^-1]` |
| `r` | Radius `[cm]` |
| `r_0` | Initial droplet radius |
| `r_s` | Instantaneous droplet radius; `r_s* = r_s/r_0` |
| `Re` | Reynolds number |
| `Sc` | Schmidt number |
| `Sh` | Sherwood number |
| `T` | Temperature `[K]` |
| `t` | Time |
| `u` | Droplet velocity |
| `u_s` | Maximum velocity on droplet surface |
| `u_∞` | Gas flow velocity |
| `v_r, v_θ` | Radial and angular liquid-velocity components |
| `W_i` | Molecular weight of component `i` |
| `X` | Distance travelled by the droplet |
| `x` | Molar fraction |
| `Y` | Mass fraction |

### Greek Symbols
| Symbol | Meaning |
|---|---|
| `α` | Thermal diffusivity |
| `δ` | Film thickness |
| `μ` | Dynamic viscosity |
| `ρ` | Density |
| `φ` | Non-dimensional parameter defined by Eq. (22) |
| `χ` | Effective conductivity parameter |

### Subscripts
| Subscript | Meaning |
|---|---|
| `air` | Air |
| `f` | Film |
| `F` | Fuel vapor |
| `g` | Gas |
| `L` | Liquid |
| `M` | Mass diffusion problem |
| `s` | Surface |
| `T` | Thermal problem |
| `∞` | Far from a droplet |
| `0` | Initial state (also without Stefan flow) |

### Superscripts / Overscore
| Mark | Meaning |
|---|---|
| `¯` | Average (reference) value |
| `*` | Modified value |
| `'` | Non-dimensional value |

## Gas-phase Transport & Film Theory

**Eq. (1) - Droplet position rate** (Revised OK)
```julia
dX_dt = u
```

**Eq. (2) - Droplet velocity rate (drag)** (Revised OK)
```julia
du_dt = (3/2*C_D/r_s*ρ_g/ρ_L)*(U_inf-u)*abs(U_inf-u)
```

**Eq. (3) - Radius regression rate** (Revised OK)
```julia
drs_dt = -ṁ/(4*π*ρ_L*r_s^2)
```

**Eq. (4) - Film reference state** (Revised OK)
```julia
T_f = T_s+A_r*(T_∞-T_s)
Y_F_f = Y_Fs+A_r*(Y_F∞-Y_Fs)
```
where A_r is the averaging parameter 

**Eq. (5) - Standard drag correlation** (Revised OK)
```julia
C_D = 24/Re*(1+Re^(2/3)/6)
```

**Eq. (6) - Thermal/diffusive film thickness** (Revised OK)
```julia
δ_T = 2*r_s/(Nu0-2)
δ_M = 2*r_s/(Sh0-2)
```

**Eq. (7) - Stefan-flow film correction factors** (Revised OK)
```julia
F_T = δ_T/δ_T0
F_M = δ_M/δ_M0
```

**Eq. (8) - Mass vaporization rate** (Revised OK)
```julia
m_dot = 2*π*ρ_g_film*r_s*D_g_film*Sh_star*log(1+B_M)
```

**Eq. (9) - Heat vaporization rate** (Revised OK)
```julia
m_dot = 2*π*k_g_film*r_s/c_p_film*Nu_star*log(1+B_T)
```
where ρ_g_film, D_g_film, k_g_film and c_p_film are the mean density, diffusion coefficient, thermal conductivity and vapor specific heat (fuel) in the film

**Eqs. (10)-(11) - Modified transfer numbers** (Revised OK)
```julia
Sh_star = 2+(Sh0-2)/F_M
Nu_star = 2+(Nu0-2)/F_T
```

**Eqs. (12)-(13) - Spalding numbers** (Revised OK)
```julia
B_M = (Y_Fs-Y_F∞)/(1-Y_Fs)
B_T = c_pF*(T_∞-T_s)/L(T_s)
```
where Y_Fs is the mass fraction of fuel vapor at the droplet surface, Y_F∞ is the mass fraction of fuel vapor far from the droplet, c_pF is the specific heat of fuel vapor, T_∞ is the temperature far from the droplet, T_s is the temperature at the droplet surface, and L(T_s) is the latent heat of vaporization at the droplet surface temperature.

**Eq. (15) - Sherwood number with blowing** (Revised OK)
```julia
Sh = Sh_star*log(1+B_M)/B_M
```

**Eqs. (16)-(17) - Universal film function** (Revised OK)
```julia
F(B) = ((1+B)^0.7*log(1+B))/B
```
Note that F(B) is a function of B, which is the Spalding number. It is used to correct the Nusselt and Sherwood numbers for the effect of blowing. and increases from 1 to 1.285 as B goes from 0 to 8. between 8 and 20 F(B) remains almost constant (1.285)
**Eq. (18) - Surface vapor composition** (Revised OK)
```julia
X_Fs = P_Fs/P
Y_Fs = X_Fs*W_F/(X_Fs*W_F+(1-X_Fs)*W_air)
```
where P_Fs is the partial pressure of fuel vapor at the droplet surface, P (revisar si es P_inf) is the partial pressure of fuel vapor far from the droplet, W_F is the molecular weight of fuel, and W_air is the molecular weight of air.

**Eq. (20) - Frossling Nu/Sh correlations** (Revised OK)
```julia
Nu0 = 2+0.552*sqrt(Re)*Pr^(1/3)
Sh0 = 2+0.552*sqrt(Re)*Sc^(1/3)
```

**Another approach Clift et al. [24]** (Revised OK)
```julia
Nu0 = 1.0+(1+Re*Pr)^(1/3)*f(Re)
Sh0 = 1.0+(1+Re*Sc)^(1/3)*f(Re)
```
where f(Re) =1 for Re<=1 and f(Re)=Re^0.077 for Re<400

**Eqs. (22)-(23) - Heat-transfer-number iteration** (Revised OK)
```julia
φ = (c_pF_film/c_p_g_film)*(Sh_star/Nu_star)*(1/Le)
B_T = (1+B_M)^(φ)-1
```

**Eq. (24) - Heat into liquid phase** (Revised OK)
```julia
Q_L = m_dot*(cp_f_film*(T_inf-T_s)/B_T+L(T_s))
```

## Internal Circulation, Reynolds & Peclet

**Eq. (25) - Liquid Reynolds and Peclet** (Revised OK)
```julia
Re_L = 2*U_s*r_s*ρ_L/μ_L
Pe_L = Re_L*Pr_L
```

**Eqs. (26)-(27) - Hill-vortex velocity field** (Revised OK)
```julia
V_r = U_s*(1-(r/r_s)^2)*cos(θ)
V_θ = U_s*(1-2*(r/r_s)^2)*sin(θ)
```

**Eq. (30) - Surface circulation speed** (Revised OK)
```julia
U_s = 1/32*(U_∞-u)*Re_inf*C_f
```
where U_s is the surface circulation speed, Re_inf is the Reynolds number based on the freestream velocity, C_f is the friction drag coefficient, U_∞ is the freestream velocity, and u is the droplet velocity.

**Eq. (31) - Friction drag with blowing** (Revised OK)
```julia
C_f = 12.69*(Re^(-2/3))
```
for 10<=Re<=100

**Eq. (32) - Friction drag with blowing in a evaporating sphere** (Revised OK)

```julia
C_f = 12.69/(Re^(2/3)*(1+B_M))
```

## Heat Equation Inside the Droplet
These are partial differential equations, not included here

## Effective Conductivity Model

**Eq. (38) - Effective conductivity definition**
```julia
χ = Nu_internal(Pe_L)/Nu_internal(0)
```

**Eq. (39) - Conductivity closure vs Peclet**
```julia
χ = 1.86+0.86*tanh(2.245*log10(Pe_L/30))
```

## Appendix — Auxiliary Mixture & Wilke Rules

### Ideal Mixture Density

**Appendix (unnumbered) - Ideal-mixture density**

```julia
ρ_mix = 1/(Y_F/ρ_F+(1-Y_F)/ρ_air)
```

### Ideal Mixture Heat Capacity
**Appendix (unnumbered) - Ideal-mixture heat capacity**
```julia
cp_mix = cp_F*Y_F+cp_air*(1-Y_F)
```

### Wilke Mixture Viscosity
**Appendix (unnumbered) - Wilke viscosity mixing**
```julia
μ_mix = sum(x_i*μ_i/(sum(x_j*Φ_ij)))
```

Where:

**Appendix (unnumbered) - Wilke interaction term**
```julia
Φ_ij = (1/√8)*(1+√(μ_i/μ_j)*(M_j/M_i)^(1/4))^2/√(1+M_i/M_j)
```

### Mixture Thermal Conductivity (Wilke)
**Appendix (unnumbered) - Wilke conductivity mixing**
```julia
k_mix = sum(x_i*k_i/(sum(x_j*Φ_ij)))
```

### Binary Diffusion (Fuller-like, simplified symbolic)
**Appendix (unnumbered) - Binary diffusion estimate**
```julia
D_F_air = D_ref*T^1.75/(P*(V_F^(1/3)+V_air^(1/3))^2)
```

---