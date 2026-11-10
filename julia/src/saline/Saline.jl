"""
    Saline

Aqueous electrolyte, saline, and water-NaCl droplet evaporation and precipitation modeling.
Implements the Raoult-van 't Hoff vapor pressure depression and salt crystallization kinetics
of Maruyama & Hasegawa (2020).
"""
module Saline

include("types.jl")
include("solver.jl")

export SalineDropletConfig, SalineSolution
export saline_vapor_pressure, saline_evaporation_constant, saline_dry_crystal_diameter
export solve_saline_droplet

end # module Saline
