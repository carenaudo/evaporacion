"""
    Types for saline droplet evaporation and salt crystallization.
    Based on Maruyama & Hasegawa (2020).
"""

"""
    SalineDropletConfig{T<:Real}

Configuration for water-salt (e.g. NaCl) droplet evaporation and precipitation.

# Fields
- `d0::T`: Initial droplet equivalent diameter [m]
- `T_sur::T`: Droplet surface temperature [K]
- `T_inf::T`: Ambient gas temperature [K]
- `RH::T`: Ambient relative humidity [%] (0 to 100)
- `P_tot::T`: Ambient total pressure [Pa] (default: 101325.0)
- `Y_s0::T`: Initial salt mass fraction [0, 1] (e.g., 0.10 for 10 wt% NaCl)
- `i_van_t_hoff::T`: Van 't Hoff dissociation factor (default: 2.0 for NaCl -> Na+ + Cl-)
- `solubility_limit::T`: Saturation salt mass fraction (default: 0.2647 for NaCl in water)
- `salt_molar_mass::T`: Molar mass of salt [kg/mol] (default: 0.05844)
- `water_molar_mass::T`: Molar mass of water [kg/mol] (default: 0.018015)
- `rho_water::T`: Pure water density [kg/m^3] (default: 997.0)
- `rho_salt::T`: Solid crystal density [kg/m^3] (default: 2160.0 for NaCl)
- `diffusivity::T`: Vapor diffusivity in air [m^2/s] (default: 2.5e-5)
"""
Base.@kwdef struct SalineDropletConfig{T<:Real}
    d0::T
    T_sur::T
    T_inf::T
    RH::T
    P_tot::T = 101325.0
    Y_s0::T = 0.10
    i_van_t_hoff::T = 2.0
    solubility_limit::T = 0.2647
    salt_molar_mass::T = 0.05844
    water_molar_mass::T = 0.018015
    rho_water::T = 997.0
    rho_salt::T = 2160.0
    diffusivity::T = 2.5e-5
end

"""
    SalineSolution{T<:Real}

Transient solution for saline droplet evaporation.
"""
struct SalineSolution{T<:Real}
    t::Vector{T}
    d::Vector{T}
    Y_s::Vector{T}
    beta::Vector{T}
    precipitation_time::T
    dry_salt_diameter::T
end
