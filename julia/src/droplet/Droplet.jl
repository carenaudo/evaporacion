"""
    Evaporacion.Droplet

Rigorous spherical droplet evaporation model with internal circulation and film blowing closures.
"""
module Droplet

using OrdinaryDiffEq
using DataFrames
using ForwardDiff

using ..Common

include("types.jl")
include("algebraic.jl")
include("surrogate_ts.jl")
include("solver.jl")
include("postprocess.jl")

export DropletParams, DropletAlgebraic, DropletLiquidState, ChebyshevTsSurrogate, BicubicLutTsSurrogate
export calc_droplet_algebraic, surface_energy_residual, intensive_surface_residual,
       precompute_liquid_state, solve_surface_temperature, solve_droplet_state,
       solve_surface_temperature_brent, solve_surface_temperature_newton,
       solve_surface_temperature_halley, solve_surface_temperature_analytical_biot,
       solve_surface_temperature_analytical_halley,
       build_chebyshev_ts_surrogate, build_lut_ts_surrogate,
       eval_chebyshev_ts, eval_lut_ts, get_or_build_ts_surrogate, clear_ts_surrogate_cache!
export droplet_ode_rhs!, droplet_dae_rhs!, make_mass_zero_callback, solve_droplet
export postprocess_trajectory
export ventilated_wet_bulb_temperature, calc_d2_revised_rate,
       reconstruct_zhou_cubic_profile, sazhin_eigenvalue_temperature_profile

end # module Droplet
