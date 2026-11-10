"""
    Evaporacion.Neural

SciML Neural ODE surrogates in Lux.jl, OrdinaryDiffEq.jl, and SciMLSensitivity.jl.
Includes Option 3b (Pure MLP Neural ODE) and Option 3c (Hybrid Regime-Attention Neural ODE).
"""
module Neural

include("option3b_mlp.jl")
include("option3c_hybrid_attn.jl")
include("checkpoint_io.jl")
include("train.jl")

using .Option3bMLP
using .Option3cHybridAttn
using .CheckpointIO
using .Train

export Option3bMLP, Option3cHybridAttn, CheckpointIO, Train
export Option3bNeuralODE, build_option3b_model, solve_neural_ode_3b
export Option3cHybridODE, build_option3c_model, solve_neural_ode_3c
export load_option3b_weights_from_npz, save_lux_weights_npz
export train_neural_ode_step!

end # module Neural
