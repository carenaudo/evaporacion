"""
    Evaporacion.Neural.Train

Native Julia training routines for Neural ODE models using Zygote and Optimisers.jl.
"""
module Train

using Lux
using Random
using Optimisers
using Zygote
using OrdinaryDiffEq
try
    using SciMLSensitivity
catch
end
using ..Option3bMLP: Option3bNeuralODE, solve_neural_ode_3b

export train_neural_ode_step!, train_epoch!

"""
    compute_loss(model, ps, st, batch_u0, batch_cond, tsteps, target_trajectories)

Computes Mean Squared Error (MSE) loss across a batch of trajectories.
"""
function compute_loss(model::Option3bNeuralODE, ps, st, batch_u0, batch_cond, tsteps, target_trajectories)
    loss = 0.0f0
    B = size(batch_u0, 2)
    for i in 1:B
        u0_i = batch_u0[:, i]
        cond_i = batch_cond[:, i]
        target_i = target_trajectories[:, :, i]

        # Integrate trajectory
        sol = solve_neural_ode_3b(model, ps, st, u0_i, cond_i, tsteps; sensealg=InterpolatingAdjoint())
        pred_i = reduce(hcat, sol.u)
        loss += sum(abs2, pred_i .- target_i) / length(target_i)
    end
    return loss / B
end

"""
    train_neural_ode_step!(model, ps, st, opt_state, batch_u0, batch_cond, tsteps, target_trajectories)

Executes a single gradient update step via Zygote reverse-mode AD and Optimisers.
"""
function train_neural_ode_step!(model::Option3bNeuralODE, ps, st, opt_state, batch_u0, batch_cond, tsteps, target_trajectories)
    loss, grads = Zygote.withgradient(p -> compute_loss(model, p, st, batch_u0, batch_cond, tsteps, target_trajectories), ps)
    opt_state, ps = Optimisers.update(opt_state, ps, grads[1])
    return loss, ps, opt_state
end

end # module Train
