"""
    Evaporacion.Neural.Option3bMLP

Lux.jl implementation of the Option 3b Pure MLP Neural ODE for [m(t), T_bar(t)].
Architecture:
    ConditionEncoder: 5 -> 64 -> 64 -> 32 (SiLU)
    NeuralRHS: (2 + 32) -> 128 -> 128 -> 64 -> 2 (SiLU)
"""
module Option3bMLP

using Lux
using Random
using OrdinaryDiffEq

export Option3bNeuralODE, build_option3b_model, solve_neural_ode_3b

"""
    build_option3b_model(; n_cond=5, emb_dim=32, hidden=128)

Builds the Lux Neural ODE model consisting of ConditionEncoder and NeuralRHS.
"""
function build_option3b_model(; n_cond::Int=5, emb_dim::Int=32, hidden::Int=128)
    cond_encoder = Chain(
        Dense(n_cond => 64, swish),
        Dense(64 => 64, swish),
        Dense(64 => emb_dim)
    )

    rhs_net = Chain(
        Dense(2 + emb_dim => hidden, swish),
        Dense(hidden => hidden, swish),
        Dense(hidden => hidden ÷ 2, swish),
        Dense(hidden ÷ 2 => 2)
    )

    return (encoder=cond_encoder, rhs=rhs_net)
end

struct Option3bNeuralODE{E, R} <: Lux.AbstractLuxContainerLayer{(:encoder, :rhs)}
    encoder::E
    rhs::R
end

function Option3bNeuralODE(; n_cond::Int=5, emb_dim::Int=32, hidden::Int=128)
    model = build_option3b_model(n_cond=n_cond, emb_dim=emb_dim, hidden=hidden)
    return Option3bNeuralODE(model.encoder, model.rhs)
end

"""
    neural_rhs_eval(model, ps, st, state_norm, cond_emb)

Evaluates the normalized derivatives [d(m_norm)/d(t_norm), d(T_norm)/d(t_norm)].
"""
function neural_rhs_eval(model::Option3bNeuralODE, ps, st, state_norm, cond_emb)
    inp = vcat(state_norm, cond_emb)
    deriv, _ = model.rhs(inp, ps.rhs, st.rhs)
    return deriv
end

"""
    solve_neural_ode_3b(model, ps, st, u0_norm, cond_norm, tsteps; solver=Tsit5(), kwargs...)

Integrates normalized state u_norm(t_norm) for t_norm ∈ [0, 1].
"""
function solve_neural_ode_3b(
    model::Option3bNeuralODE, ps, st,
    u0_norm::AbstractVector, cond_norm::AbstractVector, tsteps::AbstractVector;
    solver=Tsit5(), kwargs...
)
    # Pre-encode conditions
    cond_emb, _ = model.encoder(cond_norm, ps.encoder, st.encoder)

    f = function (u, p, t)
        inp = vcat(u, cond_emb)
        du, _ = model.rhs(inp, p, st.rhs)
        return du
    end

    prob = ODEProblem(f, u0_norm, (tsteps[1], tsteps[end]))
    sol = solve(prob, solver; p=ps.rhs, saveat=tsteps, kwargs...)
    return sol
end

end # module Option3bMLP
