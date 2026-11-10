"""
    Evaporacion.Common

Shared physical constants, High-Fidelity properties, correlations, and geometry.
"""
module Common

include("constants.jl")
include("properties.jl")
include("correlations.jl")
include("geometry.jl")

using .Constants
using .Properties
using .Correlations
using .Geometry

export Constants, Properties, Correlations, Geometry

# Re-export key functions
export R_UNIVERSAL, R_GAS, M_AIR, M_WATER, P_ATM, T_ZERO_C
export psat_water, liquid_density_water, liquid_cp_water, latent_heat_water,
       liquid_viscosity_water, liquid_thermal_conductivity_water, liquid_thermal_expansion_water,
       gas_viscosity_air, gas_thermal_conductivity_air, gas_diffusivity_water_air,
       gas_density_air, gas_cp_air, vapor_cp_water,
       vapor_viscosity_water, vapor_thermal_conductivity_water,
       wilke_mixture_viscosity, wilke_mixture_conductivity,
       reference_state_one_third, mole_fraction_from_mass_fraction, mass_fraction_from_mole_fraction,
       gas_mixture_properties,
       surface_mass_fraction_from_psat, stull_twb
export calc_re, calc_pr, calc_sc,
       nusselt_froessling, sherwood_froessling, clift_correlation,
       log1p_over_x, blowing_factor, modified_transfer_number, TransferSolution, solve_as1989_transfer,
       internal_circulation_enhancement, calc_terminal_velocity,
       spheroid_deformation_evaporation_ratio,
       holterman_evaporation_rate, holterman_relaxation_time,
       renaudo_evaporation_parameter, renaudo_extinction_diameter,
       renaudo_deposition_time, renaudo_landing_diameter,
       chiu_group_evaporation_number, chiu_evaporation_regime
export sphere_surface_from_d, sphere_volume_from_d, sphere_diameter_from_volume,
       sessile_volume_from_r, sessile_surface_from_r, r_from_volume_sessile,
       dV_dtheta, dV_dR, g_theta_exact, g_theta_fast

end # module Common
