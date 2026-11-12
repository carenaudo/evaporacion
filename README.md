# Evaporacion.jl

`Evaporacion.jl` is a Julia package containing implementations of several
droplet evaporation models. It includes models for spherical droplets in air,
sessile droplets on a surface, and some extensions for kinetic effects,
multicomponent droplets, and saline droplets. The repository also contains
scripts used to compare models and reproduce selected cases.

The implementations use different assumptions and approximations; they are
not interchangeable or necessarily complete reproductions of the cited papers.
Check the model equations, input units, and applicable regime before using
results in research or engineering calculations.

## Models

| Model family | Implementation | Main literature sources |
| --- | --- | --- |
| Spherical droplets | Coupled mass and temperature equations with Abramzon-Sirignano transfer, simpler diffusion and convection variants, and a revised D²-law variant. `DropletParams.model_variant` selects the formulation. | Abramzon & Sirignano (1989); Zhou et al. (2017); Dalla Barba et al. (2021). |
| Sessile droplets | Diffusion-driven evaporation with constant contact radius (`:CR`), constant contact angle (`:CA`), stick-slide (`:SS`), and stick-jump (`:SJ`) modes. Also includes flux, lifetime, substrate-cooling, and deposition calculations. | Wilson & D'Ambrosio (2023); Picknett & Bexon (1977); Stauber et al. (2015). |
| Kinetic effects | Diffusion-kinetic evaporation, Kelvin curvature, and corrections for interfacial and Stefan-flow effects. | Vlasov (2021); Zhao & Nadal (2023); Akkuş (2020). |
| Multicomponent droplets | Stefan-Maxwell and Stefan-Fuchs flux calculations, with solute enrichment and crust-related calculations. | Tonini & Cossali (2015); Ordoubadi et al. (2019); Gerken et al. (2014). |
| Saline droplets | Water-NaCl evaporation with vapor-pressure depression and crystallization calculations. | Maruyama & Hasegawa (2020). |
| Radial temperature model | A one-dimensional moving-boundary radial conduction model, with optional flow and surfactant effects. | See the [implementation](julia/src/radial_pde/solver.jl) for the equations and assumptions. |
| Surrogate models | Empirical volume models, surface-temperature approximations, and optional neural surrogates for comparing or approximating other models. | See [JAIIO](JAIIO/) for the repository's benchmark and validation reports. |

These are model families, not a claim that every equation or experiment in
each reference is implemented or validated. The source code and tests are the
most direct references for the available options.

## Getting started

The package lives in `julia/`. From the repository root, instantiate the Julia
environment:

```sh
julia --project=julia -e "using Pkg; Pkg.instantiate()"
```

Run the package tests with:

```sh
julia --project=julia -e "using Pkg; Pkg.test()"
```

Start a Julia session with `julia --project=julia` and try a spherical water
droplet:

```julia
using Evaporacion

params = DropletParams(
    T_inf = 35.0,                         # ambient temperature, °C
    U_rel = 0.5,                          # relative air speed, m/s
    Y_inf = 0.015,                        # ambient vapor mass fraction
    model_variant = :abramzon_sirignano
)

T0 = 22.0                                 # initial droplet temperature, °C
d0 = 1.0e-3                               # initial diameter, m
m0 = liquid_density_water(T0 + 273.15) * (π / 6) * d0^3
sol = solve_droplet(params, [m0, T0], (0.0, 400.0))
trajectory = postprocess_trajectory(sol, params)
println("Final simulated time: ", sol.t[end], " s")
println("Final mass: ", trajectory.m[end], " kg")
```

The solver can stop at its mass threshold or at the end of the requested
time span. The final time is therefore not necessarily the evaporation
lifetime.

A constant-radius sessile-droplet example:

```julia
using Evaporacion

cfg = SessileDropletConfig(
    R0 = 1.0e-3,    # initial contact radius, m
    theta0 = π / 2, # initial contact angle, radians
    mode = :CR
)

println("Analytical CR lifetime: ", lifetime_cr(cfg), " s")
```

Other inputs use model-specific units and defaults; consult the relevant
configuration type and solver in [`julia/src/`](julia/src/).

## Repository

- [`julia/`](julia/) contains the package, scripts, and tests.
- [`clicap/`](clicap/) contains scripts for selected in-flight and sessile
  droplet cases.
- [`JAIIO/`](JAIIO/) contains surrogate-model experiments and reports.
- [`docs/`](docs/) contains model notes and comparisons. The
  [models comparison](docs/MODELS_COMPARISON_AND_CHECKLIST.md) and
  [rigorous model specification](docs/RIGOROUS_MODEL_SPEC.md) provide
  additional context.

## Selected references

- Abramzon, B. and Sirignano, W. A. (1989). "Droplet vaporization model for
  spray combustion calculations." *International Journal of Heat and Mass
  Transfer*, 32, 1605-1618.
- Zhou, Z.-F., Li, W.-Y., Chen, B. and Wang, G.-X. (2017). "A 3rd-order
  polynomial temperature profile model for the heating and evaporation of
  moving droplets." *Applied Thermal Engineering*, 110, 162-170.
- Dalla Barba, F., Wang, J. and Picano, F. (2021). "Revisiting D²-law for
  the evaporation of dilute droplets." *Physics of Fluids*, 33, 051701.
- Wilson, S. K. and D'Ambrosio, H.-M. (2023). "Evaporation of sessile
  droplets." *Annual Review of Fluid Mechanics*, 55, 481-509.
- Picknett, R. G. and Bexon, R. (1977). "The evaporation of sessile or
  pendant drops in still air." *Journal of Colloid and Interface Science*,
  61, 336-350.
- Stauber, J. M., Wilson, S. K., Duffy, B. R. and Sefiane, K. (2015).
  "Evaporation of droplets on strongly hydrophobic substrates." *Langmuir*,
  31, 3653-3660.
- Vlasov, V. A. (2021). "On a theory of mass transfer during the evaporation
  of a spherical droplet." *International Journal of Heat and Mass
  Transfer*, 178, 121597.
- Zhao, H. and Nadal, F. (2023). "Droplet evaporation in inert gases."
  *Journal of Fluid Mechanics*, 958, A18.
- Akkuş, Y. (2020). "The effect of Stefan flow on the models of droplet
  evaporation." *Journal of Thermal Science and Technology*, 40, 309-318.
- Tonini, S. and Cossali, G. E. (2015). "A multi-component drop evaporation
  model based on analytical solution of Stefan-Maxwell equations."
  *International Journal of Heat and Mass Transfer*, 92, 184-189.
- Ordoubadi, M. et al. (2019). "Multi-solvent microdroplet evaporation:
  Modeling and measurement of spray-drying kinetics with inhalable
  pharmaceutics." *Pharmaceutical Research*, 36, 100.
- Gerken, W. J., Thomas, A. V., Koratkar, N. and Oehlschlaeger, M. A. (2014).
  "Nanofluid pendant droplet evaporation: Experiments and modeling."
  *International Journal of Heat and Mass Transfer*, 74, 263-268.
- Maruyama, Y. and Hasegawa, K. (2020). "Evaporation and drying kinetics of
  water-NaCl droplets via acoustic levitation." *RSC Advances*, 10, 1870-1877.
