"""
    Types for multicomponent and multi-solvent droplet evaporation.
"""

"""
    MulticomponentSpecies{T<:Real}

Physical properties of an evaporating species or non-volatile solute.

# Fields
- `name::String`: Species identifier
- `molar_mass::T`: Molar mass [kg/mol]
- `density::T`: Pure liquid density [kg/m^3]
- `boiling_temp::T`: Normal boiling point [K]
- `latent_heat::T`: Latent heat of vaporization [J/kg]
- `diffusivity::T`: Vapor diffusivity in carrier gas [m^2/s]
- `is_volatile::Bool`: Whether the species evaporates (false for solutes/nanoparticles)
- `liquid_diffusivity::T`: Internal liquid diffusion coefficient [m^2/s] (for VFL solute Péclet)
"""
Base.@kwdef struct MulticomponentSpecies{T<:Real}
    name::String
    molar_mass::T
    density::T
    boiling_temp::T = 373.15
    latent_heat::T = 2.26e6
    diffusivity::T = 2.5e-5
    is_volatile::Bool = true
    liquid_diffusivity::T = 1e-9
end

"""
    MulticomponentDropletConfig{T<:Real}

Configuration for multi-component or suspension droplet evaporation.
"""
Base.@kwdef struct MulticomponentDropletConfig{T<:Real}
    r0::T
    T_d::T
    T_inf::T
    P_tot::T = 101325.0
    species::Vector{MulticomponentSpecies{T}}
    initial_mass_fractions::Vector{T}
    ambient_mole_fractions::Vector{T} = Float64[]
    carrier_molar_mass::T = 0.02897
    carrier_diffusivity::T = 2.5e-5
end

"""
    MulticomponentSolution{T<:Real}

Transient solution for multicomponent droplet evaporation.
"""
struct MulticomponentSolution{T<:Real}
    t::Vector{T}
    r::Vector{T}
    masses::Matrix{T} # size: (nt, n_species)
    mass_fractions::Matrix{T}
    surface_enrichments::Matrix{T} # size: (nt, n_solutes)
    lifetime::T
    locking_diameter::T
end
