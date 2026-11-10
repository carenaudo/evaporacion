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
include("pinn_surrogates.jl")

using .Option3bMLP
using .Option3cHybridAttn
using .CheckpointIO
using .Train
using .PINNSurrogates

export Option3bMLP, Option3cHybridAttn, CheckpointIO, Train, PINNSurrogates
export Option3bNeuralODE, build_option3b_model, solve_neural_ode_3b
export Option3cHybridODE, build_option3c_model, solve_neural_ode_3c
export load_option3b_weights_from_npz, save_lux_weights_npz
export train_neural_ode_step!
export PINNCSurrogate, PINN2EqSurrogate, PIDeepONetSurrogate, ResMLPSurrogate,
       PINNSurrogateBundle, solve_pinn_c, solve_pinn_2eq,
       deeponet_eval_case, resmlp_eval_case,
       save_pinn_bundle, load_pinn_bundle, get_default_pinn_bundle

end # module Neural
