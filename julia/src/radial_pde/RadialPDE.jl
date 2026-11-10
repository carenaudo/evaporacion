"""
    Evaporacion.RadialPDE

1D Moving-boundary radial conduction PDE with Marangoni enhancement and surfactants.
"""
module RadialPDE

using OrdinaryDiffEq

using ..Common

include("solver.jl")

export RadialPDEConfig, RadialPDESolution, solve_radial_pde

end # module RadialPDE
