"""
    Evaporacion

Advanced Droplet Evaporation Simulation Suite in High-Performance Julia.
Fully ports and modernizes the Python package with Abramzon-Sirignano (1989) closures,
Kronig-Brink internal circulation, Wilson & D'Ambrosio (2023) sessile dynamics,
1D radial moving-boundary conduction PDE, calibrated surrogates, and SciML Neural ODEs.
"""
module Evaporacion

include("common/Common.jl")
include("droplet/Droplet.jl")
include("surrogates/Surrogates.jl")
include("sessile/Sessile.jl")
include("kinetic/Kinetic.jl")
include("multicomponent/Multicomponent.jl")
include("saline/Saline.jl")
include("radial_pde/RadialPDE.jl")
include("plotting.jl")
include("cli.jl")

using .Common
using .Droplet
using .Surrogates
using .Sessile
using .Kinetic
using .Multicomponent
using .Saline
using .RadialPDE
using .Plotting
using .CLI

export Common, Droplet, Surrogates, Sessile, Kinetic, Multicomponent, Saline, RadialPDE, Plotting, CLI

# Re-export core types and functions
export DropletParams, DropletAlgebraic, DropletLiquidState, calc_droplet_algebraic, solve_droplet, postprocess_trajectory
export ChebyshevTsSurrogate, BicubicLutTsSurrogate, build_chebyshev_ts_surrogate, build_lut_ts_surrogate,
       eval_chebyshev_ts, eval_lut_ts, get_or_build_ts_surrogate, clear_ts_surrogate_cache!
export ventilated_wet_bulb_temperature, calc_d2_revised_rate,
       reconstruct_zhou_cubic_profile, sazhin_eigenvalue_temperature_profile
export SurrogateConfig, solve_surrogate, run_parametric_sweep
export SessileDropletConfig, SessileSolution, SessileThermalParams, solve_sessile,
       cooling_effective_delta_c, estimate_interface_temperature_drop,
       lifetime_cr, lifetime_ca, lifetime_ss, lifetime_sj,
       lifetime_cr_thin, lifetime_ca_thin, lifetime_ss_thin, lifetime_sj_thin,
       lifetime_hydrophobic_sphere, stauber_critical_angle,
       lifetime_pair_cr, lifetime_pair_ss, lifetime_well_cr, lifetime_conduction_limited,
       picknett_bexon_capacitance, picknett_bexon_flux, bourges_monnier_flux, local_flux_neos, gas_conduction_index, cooperative_shielding_fluxes,
       saxton_effective_contact_angle, saxton_extinction_scaling, solve_saxton_spreading_evaporation,
       depth_averaged_radial_velocity_cca, marangoni_thermal_solutal_analysis
export RadialPDEConfig, RadialPDESolution, solve_radial_pde
export KineticDropletConfig, KineticSolution, solve_kinetic_droplet,
       kelvin_vapor_pressure_ratio, psat_kelvin, mean_thermal_speed,
       vlasov_kinetic_length, vlasov_evaporation_rate, vlasov_radius_trajectory, vlasov_lifetime,
       kinetic_peclet_number, zhao_nadal_lifetime, stefan_flow_enhancement
export MulticomponentSpecies, MulticomponentDropletConfig, MulticomponentSolution,
       stefan_maxwell_matrices, stefan_maxwell_residual, stefan_fuchs_multicomponent_fluxes,
       solute_peclet_number, vfl_surface_enrichment, crust_locking_diameter, fractal_agglomerate_max_packing,
       solve_multicomponent_droplet
export SalineDropletConfig, SalineSolution,
       saline_vapor_pressure, saline_evaporation_constant, saline_dry_crystal_diameter,
       solve_saline_droplet
export psat_water, liquid_density_water, latent_heat_water, stull_twb, g_theta_fast,
       gas_diffusivity_water_air, gas_density_air, gas_viscosity_air, gas_cp_air, gas_thermal_conductivity_air,
       liquid_cp_water, liquid_thermal_expansion_water, surface_mass_fraction_from_psat,
       gas_mixture_properties, wilke_mixture_viscosity, wilke_mixture_conductivity,
       spheroid_deformation_evaporation_ratio, holterman_evaporation_rate, holterman_relaxation_time,
       renaudo_evaporation_parameter, renaudo_extinction_diameter,
       renaudo_deposition_time, renaudo_landing_diameter, calc_terminal_velocity,
       chiu_group_evaporation_number, chiu_evaporation_regime,
       r_from_volume_sessile, sessile_volume_from_r, sphere_volume_from_d,
       R_GAS, R_UNIVERSAL, P_ATM, T_ZERO_C

# Optional Neural ODE support
try
    include("neural/Neural.jl")
    using .Neural
    export Neural, Option3bNeuralODE, Option3cHybridODE, solve_neural_ode_3b, solve_neural_ode_3c
catch
    # Lux/SciMLSensitivity not installed in the active environment
end

end # module Evaporacion
