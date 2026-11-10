"""
    Solute surface enrichment, Péclet numbers, and crust locking models.
    Based on Vehring et al. (2007), Ordoubadi et al. (2019), and Gerken et al. (2014).
"""

"""
    solute_peclet_number(kappa, D_solute)

Compute the instantaneous solute Péclet number:
```math
\\text{Pe}_j = \\frac{\\kappa}{8 D_j}
```
where ``\\kappa = -\\mathrm{d}(d^2)/\\mathrm{d}t`` is the droplet evaporation rate constant [m^2/s]
and ``D_j`` is the liquid-phase diffusion coefficient of solute ``j`` [m^2/s].
"""
function solute_peclet_number(kappa::Real, D_solute::Real)
    if D_solute <= 0.0
        return 0.0
    end
    return float(kappa) / (8.0 * float(D_solute))
end

"""
    vfl_surface_enrichment(Pe)

Calculate the steady-state solute surface enrichment ratio ``E_j = C_{\\text{surface}} / C_{\\text{mean}}``
from the Vehring-Ferron-Lechuga (VFL) model (valid for ``\\text{Pe} \\le 20``):
```math
E_j = 1 + \\frac{\\text{Pe}_j}{5} + \\frac{\\text{Pe}_j^2}{100} - \\frac{\\text{Pe}_j^3}{4000}
```
"""
function vfl_surface_enrichment(Pe::Real)
    pe = clamp(float(Pe), 0.0, 20.0)
    return 1.0 + pe / 5.0 + (pe^2) / 100.0 - (pe^3) / 4000.0
end

"""
    crust_locking_diameter(d0, phi_solute_0; porosity=0.0)

Calculate the critical locking diameter where a shell or crust forms from insoluble
nanoparticles or precipitating solute:
```math
d_{\\text{lock}} = d_0 \\left( \\frac{\\phi_{s,0}}{1 - \\epsilon_p} \\right)^{1/3}
```
where ``\\phi_{s,0}`` is initial solid volume fraction and ``\\epsilon_p`` is crust porosity.
"""
function crust_locking_diameter(d0::Real, phi_solute_0::Real; porosity::Real=0.0)
    eff_phi = clamp(float(phi_solute_0) / (1.0 - clamp(float(porosity), 0.0, 0.99)), 0.0, 1.0)
    return float(d0) * cbrt(eff_phi)
end

"""
    fractal_agglomerate_max_packing(r_np, R_g, d_f)

Calculate the maximum packing volume fraction of fractal nanoparticle agglomerates (Gerken 2014):
```math
\\phi_{\\text{max}} = \\left( \\frac{R_g}{r_{np}} \\right)^{d_f - 3} \\left( \\frac{\\pi}{6} \\right)
```
"""
function fractal_agglomerate_max_packing(r_np::Real, R_g::Real, d_f::Real)
    phi_internal = (float(R_g) / float(r_np))^(float(d_f) - 3.0)
    return phi_internal * (pi / 6.0)
end
