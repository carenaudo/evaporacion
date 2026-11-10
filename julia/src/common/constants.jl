"""
    Evaporacion.Common.Constants

Universal physical constants and species properties.
All quantities are in SI units (kg, m, s, K, Pa, J, W, mol).
"""
module Constants

export R_UNIVERSAL, R_GAS, M_AIR, M_WATER, P_ATM, T_ZERO_C, T_CRIT_WATER, T_BOIL_WATER,
       SIGMA_SB, GRAVITY, R_SPECIFIC_AIR, R_SPECIFIC_WATER

# Universal gas constant [J / (mol * K)]
const R_UNIVERSAL = 8.31446261815324
const R_GAS = R_UNIVERSAL

# Molecular weights [kg / mol]
const M_AIR   = 0.0289647    # Dry air
const M_WATER = 0.01801528   # Water (H2O)

# Standard reference conditions
const P_ATM       = 101325.0   # Standard atmospheric pressure [Pa]
const T_ZERO_C    = 273.15     # 0 °C in Kelvin [K]
const T_CRIT_WATER = 647.096   # Water critical temperature [K]
const T_BOIL_WATER = 373.15    # Water normal boiling point [K]

# Stefan-Boltzmann constant [W / (m^2 * K^4)]
const SIGMA_SB = 5.670374419e-8

# Gravitational acceleration [m / s^2]
const GRAVITY = 9.80665

# Specific gas constants [J / (kg * K)]
const R_SPECIFIC_AIR   = R_UNIVERSAL / M_AIR    # ~287.058 J/(kg*K)
const R_SPECIFIC_WATER = R_UNIVERSAL / M_WATER  # ~461.523 J/(kg*K)

end # module Constants

export R_UNIVERSAL, M_AIR, M_WATER, P_ATM, T_ZERO_C, T_CRIT_WATER, T_BOIL_WATER,
       SIGMA_SB, GRAVITY, R_SPECIFIC_AIR, R_SPECIFIC_WATER
