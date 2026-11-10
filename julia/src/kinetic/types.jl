"""
    Kinetic types for sub-micron and nanoscale droplet evaporation.
"""

"""
    KineticDropletConfig{T<:Real}

Configuration for droplet evaporation including non-equilibrium Knudsen layer kinetics,
Hertz-Knudsen phase change, Kelvin curvature shifts, and Stefan blowing.

# Fields
- `r0::T`: Initial droplet radius [m]
- `T_s::T`: Droplet surface temperature [K]
- `T_inf::T`: Ambient gas temperature [K]
- `P_tot::T`: Ambient total pressure [Pa] (default: 101325.0)
- `RH::T`: Ambient relative humidity [0, 1] (default: 0.0)
- `alpha::T`: Condensation / accommodation coefficient [0, 1] (default: 1.0)
- `molar_mass::T`: Vapor molar mass [kg/mol] (default: 0.018015 for water)
- `carrier_molar_mass::T`: Inert gas molar mass [kg/mol] (default: 0.02897 for air)
- `rho_l::T`: Liquid density [kg/m^3] (default: 997.0)
- `surface_tension::T`: Liquid-gas surface tension [N/m] (default: 0.0728)
- `use_kelvin::Bool`: Whether to include Kelvin nanoscale vapor pressure elevation (default: true)
- `use_stefan::Bool`: Whether to include Stefan convective blowing enhancement (default: true)
"""
Base.@kwdef struct KineticDropletConfig{T<:Real}
    r0::T
    T_s::T
    T_inf::T
    P_tot::T = 101325.0
    RH::T = 0.0
    alpha::T = 1.0
    molar_mass::T = 0.018015
    carrier_molar_mass::T = 0.02897
    rho_l::T = 997.0
    surface_tension::T = 0.0728
    use_kelvin::Bool = true
    use_stefan::Bool = true
end

"""
    KineticSolution{T<:Real}

Solution trajectory from kinetic droplet evaporation simulation.
"""
struct KineticSolution{T<:Real}
    t::Vector{T}
    r::Vector{T}
    mass::Vector{T}
    lifetime::T
    tau_continuum::T
    kinetic_factor::T
end
