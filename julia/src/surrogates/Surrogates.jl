"""
    Evaporacion.Surrogates

Empirical and physics-informed surrogate volume and coupled temperature models.
"""
module Surrogates

using OrdinaryDiffEq
using DataFrames

using ..Common

include("algebraic.jl")
include("solver.jl")
include("ensemble.jl")

export calc_reynolds_surrogate,
       calc_C_model4, calc_C_model5, calc_C_model6,
       calc_C_model6b, calc_C_model6b_refit,
       calc_C_model2b, calc_C_model7,
       calc_params_model8, calc_params_model8_dynamic
export SurrogateConfig, solve_surrogate, run_parametric_sweep

end # module Surrogates
