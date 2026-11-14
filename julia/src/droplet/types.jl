# Types for rigorous spherical droplet evaporation

export DropletParams, DropletAlgebraic, DropletLiquidState, ChebyshevTsSurrogate, BicubicLutTsSurrogate, TsSolverCounter, reset!

"""
    TsSolverCounter

Lightweight mutable counter for tracking root-finder calls and inner iterations during benchmarking.
"""
mutable struct TsSolverCounter
    calls::Int
    iters::Int
end
TsSolverCounter() = TsSolverCounter(0, 0)
reset!(c::TsSolverCounter) = (c.calls = 0; c.iters = 0; c)

"""
    ChebyshevTsSurrogate

Compact 2D (or 1D when Re=0) Chebyshev polynomial surrogate for equilibrium surface temperature Ts(rs, T_bar).
Domain: xi = sqrt(rs) in [0, xi_max], T_bar in [T_min, T_max].
"""
struct ChebyshevTsSurrogate
    xi_max::Float64
    T_min::Float64
    T_max::Float64
    is_1d::Bool
    coeffs::Matrix{Float64}  # size (n_xi, n_T) or (1, n_T)
end

"""
    BicubicLutTsSurrogate

2D Lookup Table with bicubic Hermite (Catmull-Rom) interpolation for Ts(rs, T_bar).
Domain: xi = sqrt(rs) in [0, xi_max], T_bar in [T_min, T_max].
"""
struct BicubicLutTsSurrogate
    xi_max::Float64
    T_min::Float64
    T_max::Float64
    inv_dxi::Float64
    inv_dT::Float64
    data::Matrix{Float64}  # size (n_xi, n_T)
end

"""
    DropletParams{T<:Real}

Configuration parameters for rigorous spherical droplet evaporation.
"""
Base.@kwdef struct DropletParams{T<:Real}
    T_inf::T
    U_rel::Float64 = 0.0
    Y_inf::T
    P_atm::Float64 = 101325.0
    profile::Symbol = :parabolico
    correlation::Symbol = :frossling
    frossling_coef::Float64 = 0.552             # Frössling coefficient in Sh0/Nu0 = 2 + c Re^1/2 (Sc|Pr)^1/3, used by every model variant
    as1989_blowing::Bool = true
    fixed_chi_m::Bool = false
    model_variant::Symbol = :abramzon_sirignano  # :abramzon_sirignano, :level_1, :level_2, :d2_classic, :d2_stefan, :d2_convective, :d2_revised_dallabarba
    terminal_velocity::Bool = false             # dynamically evaluate U_rel = calc_terminal_velocity(d, T_eg_K, P_atm)
    isothermal::Bool = false                    # if true, droplet temperature is held fixed (dT/dt = 0)
    T_d_fixed::Float64 = NaN                    # fixed droplet temperature (°C); if NaN, defaults to Stull T_wb
    thermal_dilation::Bool = false              # Angelilli et al. (2024) liquid thermal dilation
    variable_cp_factor::Bool = false            # Angelilli et al. (2024) variable-cpl enthalpy correction factor
    spheroid_epsilon::Float64 = 0.0             # Tonini & Cossali (2015) deformation parameter ε = a/b - 1
    zhou_moving_boundary::Bool = true           # Zhou et al. (2017) moving boundary (Ts - T_bar) term in dT/dt
    reference_state::Symbol = :average          # :average (T_eg = (Ts+Tinf)/2) or :one_third (Hubbard/A&S 1/3 rule)
    use_wilke_mixtures::Bool = false            # when true, uses Wilke binary mixing for gas viscosity and conductivity
    ts_solver::Symbol = :newton                 # :brent, :newton, :halley, :analytical_biot, :analytical_halley, :dae, :chebyshev, :lut
    ts_surrogate::Any = nothing                 # precomputed ChebyshevTsSurrogate or BicubicLutTsSurrogate (or nothing for auto-cache)
    ts_counter::Union{Nothing, TsSolverCounter} = nothing # optional counter for Ts solver calls & iterations
end

"""
    DropletLiquidState{T<:Real}

Encapsulates state properties of the droplet that depend strictly on (m, T_bar)
and are invariant to the interface temperature Ts.
"""
struct DropletLiquidState{T<:Real}
    m::T
    T_bar::T
    T_bar_K::T
    rho_l::T
    cp_l::T
    k_l::T
    mu_l::T
    d::T
    r_s::T
    u_rel::T
end

"""
    DropletAlgebraic{T<:Real}

Complete diagnostic output of the algebraic block at a given state (m, T_bar, Ts).
"""
struct DropletAlgebraic{T<:Real}
    d::T            # Diameter [m]
    r_s::T          # Droplet radius [m]
    T_eg::T         # Gas film temperature [°C]
    rho_g::T        # Gas density at film [kg/m^3]
    mu_g::T         # Gas dynamic viscosity [Pa*s]
    cp_g::T         # Gas specific heat [J/(kg*K)]
    k_g::T          # Gas thermal conductivity [W/(m*K)]
    D_ab::T         # Binary diffusivity [m^2/s]
    cp_l::T         # Liquid specific heat [J/(kg*K)]
    k_l::T          # Liquid thermal conductivity [W/(m*K)]
    rho_l::T        # Liquid density [kg/m^3]
    DHev::T         # Latent heat of vaporization [J/kg]
    Re::T           # Reynolds number [-]
    Pr::T           # Prandtl number [-]
    Sc::T           # Schmidt number [-]
    Nu0::T          # Base Nusselt number [-]
    Sh0::T          # Base Sherwood number [-]
    Y_s::T          # Surface vapor mass fraction [-]
    B_M::T          # Spalding mass number [-]
    B_T::T          # Spalding heat number [-]
    F_M::T          # Mass blowing factor [-]
    F_T::T          # Heat blowing factor [-]
    Sh_star::T      # Modified Sherwood number [-]
    Sh::T           # Effective Sherwood number [-]
    Nu_star::T      # Modified Nusselt number [-]
    Nu::T           # Effective Nusselt number [-]
    m_dot::T        # Evaporation mass rate [kg/s]
    h::T            # Convective heat transfer coefficient [W/(m^2*K)]
    Le::T           # Lewis number [-]
    r_dot_s::T      # Radius shrinkage rate [m/s]
    q_l::T          # Liquid heat conduction rate [W]
    q_g::T          # Gas convective heat transfer rate [W]
    b_q::T          # Thermal profile coefficient [-]
    b_m::T          # Mass profile coefficient [-]
    chi_l::T        # Thermal circulation enhancement factor [-]
    chi_m::T        # Mass circulation enhancement factor [-]
    U_s::T          # Surface circulation velocity [m/s]
end
