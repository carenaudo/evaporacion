"""
    Evaporacion.Neural.Option3cHybridAttn

Lux.jl implementation of the Option 3c Hybrid Attention Neural ODE.
Combines a differentiable physics baseline with a Regime-Attention transformer block
that learns residual dynamics across physical evaporation regimes.
"""
module Option3cHybridAttn

using Lux
using Random
using OrdinaryDiffEq
using LinearAlgebra
using ...Evaporacion.Droplet: DropletParams, calc_droplet_algebraic, solve_surface_temperature

export Option3cHybridODE, build_option3c_model, solve_neural_ode_3c

"""
    RegimeAttentionLayer

Cross-attention layer that attends over learnable regime prototype vectors.
"""
struct RegimeAttentionLayer <: Lux.AbstractLuxLayer
    d_model::Int
    n_regimes::Int
end

function RegimeAttentionLayer(; d_model::Int=64, n_regimes::Int=8)
    return RegimeAttentionLayer(d_model, n_regimes)
end

function Lux.initialparameters(rng::AbstractRNG, l::RegimeAttentionLayer)
    # Prototypes matrix: (d_model, n_regimes)
    prototypes = randn(rng, Float32, l.d_model, l.n_regimes) .* 0.02f0
    return (prototypes=prototypes,)
end

Lux.initialstates(::AbstractRNG, ::RegimeAttentionLayer) = (;)

function ((l::RegimeAttentionLayer)(q::AbstractVector, ps, st))
    # q is (d_model,)
    # ps.prototypes is (d_model, n_regimes)
    scale = 1.0f0 / sqrt(Float32(l.d_model))
    # Dot-product attention scores
    scores = (transpose(ps.prototypes) * q) .* scale # (n_regimes,)
    # Softmax
    scores_max = maximum(scores)
    exp_scores = exp.(scores .- scores_max)
    attn_weights = exp_scores ./ sum(exp_scores)
    # Attention context vector
    context = ps.prototypes * attn_weights # (d_model,)
    return context, st
end

"""
    build_option3c_model(; d_model=64, n_regimes=8)

Builds the Hybrid Attention network:
Query projection -> Regime cross-attention -> FeedForward -> Residual head.
"""
function build_option3c_model(; d_model::Int=64, n_regimes::Int=8)
    query_proj = Chain(
        Dense(2 + 5 => d_model * 2, gelu),
        Dense(d_model * 2 => d_model)
    )

    attn_layer = RegimeAttentionLayer(d_model=d_model, n_regimes=n_regimes)

    ff_block = Chain(
        Dense(d_model => d_model * 4, gelu),
        Dense(d_model * 4 => d_model)
    )

    residual_head = Chain(
        Dense(d_model => d_model, gelu),
        Dense(d_model => 2)
    )

    return (
        query=query_proj,
        attn=attn_layer,
        ff=ff_block,
        residual=residual_head
    )
end

struct Option3cHybridODE{Q, A, F, R} <: Lux.AbstractLuxContainerLayer{(:query, :attn, :ff, :residual)}
    query::Q
    attn::A
    ff::F
    residual::R
end

function Option3cHybridODE(; d_model::Int=64, n_regimes::Int=8)
    m = build_option3c_model(d_model=d_model, n_regimes=n_regimes)
    return Option3cHybridODE(m.query, m.attn, m.ff, m.residual)
end

"""
    evaluate_hybrid_residual(model, ps, st, state_cond)

Evaluates the attention-weighted residual [Δdm/dt, ΔdT/dt].
"""
function evaluate_hybrid_residual(model::Option3cHybridODE, ps, st, state_cond::AbstractVector)
    q, _ = model.query(state_cond, ps.query, st.query)
    attn_ctx, _ = model.attn(q, ps.attn, st.attn)
    # Residual connection + FeedForward
    x = q .+ attn_ctx
    ff_out, _ = model.ff(x, ps.ff, st.ff)
    x2 = x .+ ff_out
    delta_deriv, _ = model.residual(x2, ps.residual, st.residual)
    return delta_deriv
end

"""
    solve_neural_ode_3c(model, ps, st, phys_params, u0, cond_vec, tspan; solver=Tsit5(), kwargs...)

Integrates the hybrid physical + attention residual ODE system:
    du/dt = du_physics/dt + du_residual/dt
"""
function solve_neural_ode_3c(
    model::Option3cHybridODE, ps, st,
    phys_params::DropletParams,
    u0::AbstractVector, cond_vec::AbstractVector, tspan::Tuple{Real, Real};
    solver=Tsit5(), kwargs...
)
    f = function (u, p, t)
        m = u[1]
        T_bar = u[2]

        if m <= 1e-18
            return zeros(Float32, 2)
        end

        # 1. Physics baseline
        Ts = solve_surface_temperature(m, T_bar, phys_params)
        alg = calc_droplet_algebraic(m, T_bar, Ts, phys_params)

        dm_phys = -Float32(alg.m_dot)
        term1 = 3.0f0 * Float32(alg.h) * (Float32(phys_params.T_inf) - Float32(Ts)) / (Float32(alg.rho_l) * Float32(alg.cp_l) * Float32(alg.r_s))
        term2 = 3.0f0 * Float32(alg.r_dot_s) * Float32(alg.DHev) / (Float32(alg.r_s) * Float32(alg.cp_l))
        dT_phys = term1 + term2

        # 2. Attention residual
        state_cond = vcat(Float32[m, T_bar], Float32.(cond_vec))
        delta_u = evaluate_hybrid_residual(model, ps, st, state_cond)

        return [dm_phys + delta_u[1], dT_phys + delta_u[2]]
    end

    prob = ODEProblem(f, Float32.(u0), Float32.(tspan))
    return solve(prob, solver; kwargs...)
end

end # module Option3cHybridAttn
