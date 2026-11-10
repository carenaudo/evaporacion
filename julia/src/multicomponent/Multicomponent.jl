"""
    Multicomponent

Competitive multicomponent, multi-solvent, and suspension droplet evaporation modeling.
Implements the exact analytical Stefan-Maxwell matrix exponential model of Tonini & Cossali (2015),
the Stefan-Fuchs radial advection model of Ordoubadi et al. (2019), and Vehring-Ferron-Lechuga (VFL)
solute surface enrichment and crust locking of Gerken et al. (2014).
"""
module Multicomponent

include("types.jl")
include("stefan_maxwell.jl")
include("surface_enrichment.jl")
include("solver.jl")

export MulticomponentSpecies, MulticomponentDropletConfig, MulticomponentSolution
export stefan_maxwell_matrices, stefan_maxwell_residual, stefan_fuchs_multicomponent_fluxes
export solute_peclet_number, vfl_surface_enrichment, crust_locking_diameter, fractal_agglomerate_max_packing
export solve_multicomponent_droplet

end # module Multicomponent
