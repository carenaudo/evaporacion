"""
    Kinetic

Sub-micron, nanoscale, and Knudsen-layer non-equilibrium droplet evaporation modeling.
Implements the diffusion-kinetic theory of Vlasov (2021), non-equilibrium Knudsen
jump conditions of Zhao & Nadal (2023), and Stefan interfacial blowing of Akkuş (2020).
"""
module Kinetic

include("types.jl")
include("solver.jl")

export KineticDropletConfig, KineticSolution
export kelvin_vapor_pressure_ratio, psat_kelvin, mean_thermal_speed
export vlasov_kinetic_length, vlasov_evaporation_rate, vlasov_radius_trajectory, vlasov_lifetime
export kinetic_peclet_number, zhao_nadal_lifetime
export stefan_flow_enhancement, solve_kinetic_droplet

end # module Kinetic
