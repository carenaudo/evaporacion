"""
    evaporacion_gota.jl

Top-level runner script for the Evaporacion.jl rigorous droplet evaporation model.
Abramzon & Sirignano (1989) + Zhou et al. (2017) High-Fidelity formulation.
"""

# Activate the local package environment
import Pkg
if !haskey(Pkg.project().dependencies, "Evaporacion")
    include("src/Evaporacion.jl")
    using .Evaporacion
else
    using Evaporacion
end

# Delegate to CLI main
Evaporacion.CLI.main(ARGS)
