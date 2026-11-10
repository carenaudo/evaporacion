"""
    Evaporacion.Sessile

Sessile droplet evaporation suite for CR, CA, SS, and SJ modes,
conductive substrate cooling, internal flows, and coffee-ring deposition.
Wilson & D'Ambrosio (2023) and Popov (2005).
"""
module Sessile

using OrdinaryDiffEq
using QuadGK
using LinearAlgebra

using ..Common

include("flux.jl")
include("solver.jl")
include("internal_deposition.jl")

export SessileThermalParams, popov_global_flux, cooling_effective_delta_c,
       estimate_interface_temperature_drop, local_flux_popov_series, local_flux_deegan,
       picknett_bexon_capacitance, picknett_bexon_flux, bourges_monnier_flux,
       local_flux_neos, gas_conduction_index, cooperative_shielding_fluxes
export SessileDropletConfig, SessileSolution,
       lifetime_cr, lifetime_ca, lifetime_ss, lifetime_sj,
       lifetime_cr_thin, lifetime_ca_thin, lifetime_ss_thin, lifetime_sj_thin,
       lifetime_hydrophobic_sphere, stauber_critical_angle,
       lifetime_pair_cr, lifetime_pair_ss, lifetime_well_cr,
       lifetime_conduction_limited,
       saxton_effective_contact_angle, saxton_extinction_scaling, solve_saxton_spreading_evaporation,
       solve_sessile
export depth_averaged_radial_velocity, depth_averaged_radial_velocity_cca,
       velocity_field_u, velocity_field_w,
       coffee_ring_deposit_scaling, marangoni_thermal_solutal_analysis

end # module Sessile
