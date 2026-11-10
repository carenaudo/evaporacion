"""
    jaiio_common.jl

Shared configuration, model registry, Muller (2016) CWU complexity weights,
and unified trajectory runner for all JAIIO 2026 Julia scripts.
"""

module JAIIOCommon

using Printf
using Statistics
using LinearAlgebra
using DataFrames
using CSV
using OrdinaryDiffEq

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const JAIIO_DIR = @__DIR__
const DATA_DIR  = joinpath(JAIIO_DIR, "data")
const RESULTS_DIR = joinpath(JAIIO_DIR, "results")
const TABLES_DIR  = joinpath(RESULTS_DIR, "tables")
const FIGURES_BENCH_DIR = joinpath(JAIIO_DIR, "figures", "benchmark")
const FIGURES_VAL_DIR   = joinpath(JAIIO_DIR, "figures", "validation")
const FIGURES_PINN_DIR  = joinpath(JAIIO_DIR, "figures", "pinn")
const PINN_WEIGHTS_PATH = joinpath(RESULTS_DIR, "pinn_weights.npz")

for d in (DATA_DIR, RESULTS_DIR, TABLES_DIR, FIGURES_BENCH_DIR, FIGURES_VAL_DIR, FIGURES_PINN_DIR)
    mkpath(d)
end

const JULIA_PKG_DIR = joinpath(REPO_ROOT, "julia")
if !(JULIA_PKG_DIR in LOAD_PATH)
    push!(LOAD_PATH, JULIA_PKG_DIR)
end

using Evaporacion
using Evaporacion.Common
using Evaporacion.Droplet
using Evaporacion.Surrogates
using Evaporacion.Neural
using Evaporacion.Neural.PINNSurrogates

export REPO_ROOT, JAIIO_DIR, DATA_DIR, RESULTS_DIR, TABLES_DIR,
       FIGURES_BENCH_DIR, FIGURES_VAL_DIR, FIGURES_PINN_DIR, PINN_WEIGHTS_PATH,
       MODEL_SPECS, PINN_MODEL_SPECS, ALL_MODEL_SPECS, MODEL_MAP,
       TS_SOLVER_SPECS, MULLER_CWU_WEIGHTS,
       V2_GRID, V1_GRID, BENCHMARK_CASES_25,
       get_v2_conditions, get_v1_conditions,
       display_name, model_color, model_category, cwu_weight,
       load_surrogate_coeffs, save_surrogate_coeffs,
       load_pinn_models,
       run_jaiio_case, solve_jaiio_trajectory

# -----------------------------------------------------------------------------
# Model Metadata Registry (Consistent JAIIO Naming & Palette)
# -----------------------------------------------------------------------------

struct ModelSpec
    id::Symbol
    display::String
    category::Symbol   # :physical, :surrogate, or :neural
    color::String
    description::String
end

const MODEL_SPECS = [
    ModelSpec(:rigorous,         "Riguroso",           :physical,  "#2c3e50", "Nivel 3 A&S + perfil parabólico + solver de Ts"),
    ModelSpec(:cubic,            "Cúbico",             :physical,  "#1abc9c", "Nivel 3 A&S + perfil cúbico Zhou (2017)"),
    ModelSpec(:level_2,          "No-Stefan",          :physical,  "#34495e", "Físico reducido lineal con convección (Sh0, Nu0)"),
    ModelSpec(:level_1,          "Difusión-Pura",      :physical,  "#8e44ad", "Físico reducido difusión pura lineal (Sh=2, Nu=2)"),
    ModelSpec(:uniform,          "Perfil-Uniforme",    :physical,  "#7f8c8d", "Nivel 3 A&S con perfil uniforme (Ts = T_bar, sin root-finder)"),
    ModelSpec(:iso_twb,          "Equilibrio-Térmico", :physical,  "#95a5a6", "Gota isotérmica fija a T_wb (sin EDO térmica)"),
    ModelSpec(:modelo_2b,        "Aproximado-1",       :surrogate, "#e74c3c", "Ley d² optimizada evaluada en T_wb (ex modelo_2b)"),
    ModelSpec(:modelo_6b,        "Modelo-6b",          :surrogate, "#f39c12", "Surrogado empírico segmentado coeficientes originales"),
    ModelSpec(:modelo_6b_refit,  "Aproximado-2",       :surrogate, "#27ae60", "Surrogado empírico segmentado recalibrado (ex modelo_6b_refit)"),
    ModelSpec(:modelo_7,         "Aproximado-3",       :surrogate, "#2980b9", "Surrogado multi-variable HR y Θ recalibrado (ex modelo_7)"),
]

const PINN_MODEL_SPECS = [
    ModelSpec(:pinn_c,      "PINN-1Eq",    :neural, "#16a085", "PINN 1-EDO ley d² residual suave C∞ + penalización física de monotonicidad"),
    ModelSpec(:pinn_2eq,    "PINN-2Eq",    :neural, "#8e44ad", "PINN 2-EDO acoplada [V(t), T̄(t)] + cierre neuronal Biot Ts"),
    ModelSpec(:pi_deeponet, "PI-DeepONet", :neural, "#2c82c9", "Operador DeepONet Physics-Informed sin integrador (nfev=1)"),
    ModelSpec(:resmlp,      "ResMLP-Puro", :neural, "#7f8c8d", "Red MLP pura sin física previa (línea base ablación)"),
]

const ALL_MODEL_SPECS = vcat(MODEL_SPECS, PINN_MODEL_SPECS)

const MODEL_MAP = Dict(s.id => s for s in ALL_MODEL_SPECS)

display_name(id::Symbol) = haskey(MODEL_MAP, id) ? MODEL_MAP[id].display : String(id)
display_name(id::AbstractString) = display_name(Symbol(id))
model_color(id::Symbol) = haskey(MODEL_MAP, id) ? MODEL_MAP[id].color : "#555555"
model_color(id::AbstractString) = model_color(Symbol(id))
model_category(id::Symbol) = haskey(MODEL_MAP, id) ? MODEL_MAP[id].category : :surrogate

# New Julia Ts solvers comparison list
const TS_SOLVER_SPECS = [
    (id=:brent,             display="Brent (Ref. 1973)"),
    (id=:newton,            display="Newton-Raphson (Warm)"),
    (id=:halley,            display="Halley 3er Orden"),
    (id=:analytical_biot,   display="Biot Analítico (0-iter)"),
    (id=:analytical_halley, display="Halley Analítico (0-iter)"),
    (id=:chebyshev,         display="Chebyshev 2D Surrogate"),
    (id=:lut,               display="Bicubic Hermite LUT"),
    (id=:dae,               display="Index-1 DAE (Rodas5P)"),
]

# -----------------------------------------------------------------------------
# Muller (2016) Complexity Weights (CWU)
# -----------------------------------------------------------------------------
# Weights per RHS / root residual evaluation relative to Aproximado-2 (modelo_6b_refit)
# Baseline weights: +,*=1, /=4, sqrt=4, cbrt=6, exp/log/pow=15.
function compute_muller_cwu_weights()
    cost_phys = 185.0
    cost_m6b = 87.0
    cost_m2b = 150.0
    cost_m7  = 122.0

    raw = Dict{Symbol, Float64}(
        :rigorous => cost_phys,
        :cubic    => cost_phys,
        :level_2  => cost_phys * 0.90,
        :level_1  => cost_phys * 0.85,
        :uniform  => cost_phys,
        :iso_twb  => cost_phys * 0.95,
        :modelo_2b        => cost_m2b,
        :modelo_6b        => cost_m6b,
        :modelo_6b_refit  => cost_m6b,
        :modelo_7         => cost_m7,
        # Neural & PINN weights (transparently counting MACs + nonlinear activations):
        :pinn_c           => 1735.0,
        :pinn_2eq         => 2840.0,
        :pi_deeponet      => 14500.0,
        :resmlp           => 18000.0,
    )
    return Dict(k => round(v / cost_m6b; digits=3) for (k, v) in raw)
end

const MULLER_CWU_WEIGHTS = compute_muller_cwu_weights()
cwu_weight(id::Symbol) = get(MULLER_CWU_WEIGHTS, id, 1.0)
cwu_weight(id::AbstractString) = cwu_weight(Symbol(id))

# -----------------------------------------------------------------------------
# Parametric Study Grids & 25-Case Benchmark Set
# -----------------------------------------------------------------------------

const V2_GRID = (
    d0    = [100e-6, 200e-6, 500e-6, 1000e-6, 1500e-6, 2000e-6, 2800e-6],
    U     = [0.0, 0.5, 1.0, 2.0, 3.0],
    T0    = [10.0, 15.0, 22.0, 35.0, 50.0],
    T_inf = [15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0],
    HR    = [0.3, 0.5, 0.7],
) # 7 * 5 * 24 (T_inf > T0) * 3 = 2,520 cases

const V1_GRID = (
    d0    = [100e-6, 500e-6, 1000e-6, 2000e-6, 2800e-6],
    U     = [0.0, 0.5, 1.5, 3.0],
    T0    = [10.0, 22.0, 50.0],
    T_inf = [20.0, 35.0, 50.0, 65.0],
    HR    = [0.1, 0.4, 0.7, 0.9],
) # 5 * 4 * 8 (T_inf > T0) * 4 = 640 cases

function _build_grid_conditions(grid)
    conds = NamedTuple{(:d0_m, :U_ms, :T0_C, :T_inf_C, :HR), NTuple{5, Float64}}[]
    for d0 in grid.d0, U in grid.U, T0 in grid.T0, Tinf in grid.T_inf, HR in grid.HR
        if Tinf > T0
            push!(conds, (d0_m=d0, U_ms=U, T0_C=T0, T_inf_C=Tinf, HR=HR))
        end
    end
    return conds
end

get_v2_conditions() = _build_grid_conditions(V2_GRID)
get_v1_conditions() = _build_grid_conditions(V1_GRID)

# Exact 25 Representative Cases (4 extreme corners + 6 edge + 15 interior)
const BENCHMARK_CASES_25 = [
    # 4 extreme corners
    (case_id=1,  d0_m=100e-6,  U_ms=3.0, T0_C=10.0, T_inf_C=50.0, HR=0.3, note="smallest, fast air, hot, dry"),
    (case_id=2,  d0_m=2800e-6, U_ms=0.0, T0_C=35.0, T_inf_C=40.0, HR=0.7, note="largest, no air, humid"),
    (case_id=3,  d0_m=100e-6,  U_ms=0.0, T0_C=10.0, T_inf_C=15.0, HR=0.7, note="ultra-small, cold, humid"),
    (case_id=4,  d0_m=2800e-6, U_ms=3.0, T0_C=35.0, T_inf_C=50.0, HR=0.3, note="large, hot, convective, dry"),
    # 6 edge cases
    (case_id=5,  d0_m=100e-6,  U_ms=1.0, T0_C=22.0, T_inf_C=35.0, HR=0.5, note="smallest d0"),
    (case_id=6,  d0_m=2800e-6, U_ms=1.0, T0_C=22.0, T_inf_C=35.0, HR=0.5, note="largest d0"),
    (case_id=7,  d0_m=500e-6,  U_ms=0.0, T0_C=22.0, T_inf_C=35.0, HR=0.5, note="zero airflow"),
    (case_id=8,  d0_m=500e-6,  U_ms=3.0, T0_C=22.0, T_inf_C=35.0, HR=0.5, note="max airflow"),
    (case_id=9,  d0_m=500e-6,  U_ms=1.0, T0_C=22.0, T_inf_C=35.0, HR=0.3, note="driest"),
    (case_id=10, d0_m=500e-6,  U_ms=1.0, T0_C=22.0, T_inf_C=35.0, HR=0.7, note="most humid"),
    # 15 interior spread
    (case_id=11, d0_m=200e-6,  U_ms=0.5, T0_C=15.0, T_inf_C=25.0, HR=0.3, note="interior 1"),
    (case_id=12, d0_m=200e-6,  U_ms=2.0, T0_C=22.0, T_inf_C=40.0, HR=0.5, note="interior 2"),
    (case_id=13, d0_m=500e-6,  U_ms=1.0, T0_C=15.0, T_inf_C=30.0, HR=0.5, note="interior 3"),
    (case_id=14, d0_m=500e-6,  U_ms=2.0, T0_C=35.0, T_inf_C=45.0, HR=0.3, note="interior 4"),
    (case_id=15, d0_m=500e-6,  U_ms=0.5, T0_C=10.0, T_inf_C=25.0, HR=0.7, note="interior 5"),
    (case_id=16, d0_m=1000e-6, U_ms=0.0, T0_C=22.0, T_inf_C=30.0, HR=0.5, note="interior 6"),
    (case_id=17, d0_m=1000e-6, U_ms=1.0, T0_C=15.0, T_inf_C=40.0, HR=0.3, note="interior 7"),
    (case_id=18, d0_m=1000e-6, U_ms=2.0, T0_C=35.0, T_inf_C=50.0, HR=0.7, note="interior 8"),
    (case_id=19, d0_m=1000e-6, U_ms=3.0, T0_C=22.0, T_inf_C=45.0, HR=0.5, note="interior 9"),
    (case_id=20, d0_m=1500e-6, U_ms=0.5, T0_C=22.0, T_inf_C=35.0, HR=0.3, note="interior 10"),
    (case_id=21, d0_m=1500e-6, U_ms=2.0, T0_C=15.0, T_inf_C=30.0, HR=0.7, note="interior 11"),
    (case_id=22, d0_m=2000e-6, U_ms=1.0, T0_C=35.0, T_inf_C=45.0, HR=0.5, note="interior 12"),
    (case_id=23, d0_m=2000e-6, U_ms=0.5, T0_C=10.0, T_inf_C=20.0, HR=0.3, note="interior 13"),
    (case_id=24, d0_m=2000e-6, U_ms=3.0, T0_C=22.0, T_inf_C=40.0, HR=0.7, note="interior 14"),
    (case_id=25, d0_m=200e-6,  U_ms=3.0, T0_C=50.0, T_inf_C=50.0, HR=0.3, note="warm start T0=Tinf"),
]

# -----------------------------------------------------------------------------
# Fitted Surrogate & PINN Coefficients Persistence
# -----------------------------------------------------------------------------

const FITTED_COEFFS_PATH = joinpath(RESULTS_DIR, "fitted_surrogate_coeffs.csv")

function save_surrogate_coeffs(m6b_coeffs, m7_coeffs; path::AbstractString=FITTED_COEFFS_PATH)
    rows = NamedTuple[]
    for (reg, vals) in ((:low, m6b_coeffs.low), (:mid, m6b_coeffs.mid), (:high, m6b_coeffs.high))
        push!(rows, (model="modelo_6b_refit", regime=string(reg),
                     A0=vals[1], n0=vals[2], alpha=0.0,
                     A1=vals[3], n1=vals[4], beta=0.0,
                     gamma=0.0, delta=0.0, ci=vals[5]))
    end
    for (reg, vals) in ((:low, m7_coeffs.low), (:mid, m7_coeffs.mid), (:high, m7_coeffs.high))
        push!(rows, (model="modelo_7", regime=string(reg),
                     A0=vals[1], n0=vals[2], alpha=vals[3],
                     A1=vals[4], n1=vals[5], beta=vals[6],
                     gamma=vals[7], delta=vals[8], ci=vals[9]))
    end
    df = DataFrame(rows)
    CSV.write(path, df)
    return path
end

function load_surrogate_coeffs(; path::AbstractString=FITTED_COEFFS_PATH)
    if !isfile(path)
        return (modelo_6b_refit=nothing, modelo_7=nothing)
    end
    df = CSV.read(path, DataFrame)

    function extract_m6b(sub)
        d = Dict(Symbol(r.regime) => (r.A0, r.n0, r.A1, r.n1, r.ci) for r in eachrow(sub))
        (low=d[:low], mid=d[:mid], high=d[:high])
    end

    function extract_m7(sub)
        d = Dict(Symbol(r.regime) => (r.A0, r.n0, r.alpha, r.A1, r.n1, r.beta, r.gamma, r.delta, r.ci) for r in eachrow(sub))
        (low=d[:low], mid=d[:mid], high=d[:high])
    end

    sub6 = filter(r -> r.model == "modelo_6b_refit", df)
    sub7 = filter(r -> r.model == "modelo_7", df)
    return (
        modelo_6b_refit = nrow(sub6) == 3 ? extract_m6b(sub6) : nothing,
        modelo_7        = nrow(sub7) == 3 ? extract_m7(sub7)  : nothing,
    )
end

function load_pinn_models(; path::AbstractString=PINN_WEIGHTS_PATH)
    return get_default_pinn_bundle(path)
end

# -----------------------------------------------------------------------------
# Unified Case Runner (Physical + Surrogate + PINN/Neural)
# -----------------------------------------------------------------------------

@inline function _extract_nfev(sol)
    if hasproperty(sol, :stats)
        return Int(sol.stats.nf)
    elseif hasproperty(sol, :destats)
        return Int(sol.destats.nf)
    elseif hasproperty(sol, :nfev)
        return Int(sol.nfev)
    end
    return length(sol.t)
end

function build_physical_params(
    model_id::Symbol, d0_m::Float64, U_ms::Float64, T0_C::Float64, T_inf_C::Float64, HR::Float64;
    ts_solver::Symbol=:brent, ts_counter::Union{Nothing, TsSolverCounter}=nothing
)
    psat_inf = psat_water(T_inf_C + T_ZERO_C)
    y_inf    = surface_mass_fraction_from_psat(HR * psat_inf, P_ATM)
    twb      = stull_twb(T_inf_C, HR * 100.0)

    if model_id === :rigorous
        return DropletParams(T_inf=T_inf_C, U_rel=U_ms, Y_inf=y_inf, P_atm=P_ATM,
                             profile=:parabolico, as1989_blowing=true, model_variant=:abramzon_sirignano,
                             ts_solver=ts_solver, ts_counter=ts_counter)
    elseif model_id === :cubic
        return DropletParams(T_inf=T_inf_C, U_rel=U_ms, Y_inf=y_inf, P_atm=P_ATM,
                             profile=:cubico, as1989_blowing=true, model_variant=:abramzon_sirignano,
                             ts_solver=ts_solver, ts_counter=ts_counter)
    elseif model_id === :level_2
        return DropletParams(T_inf=T_inf_C, U_rel=U_ms, Y_inf=y_inf, P_atm=P_ATM,
                             profile=:parabolico, as1989_blowing=false, model_variant=:level_2,
                             ts_solver=ts_solver, ts_counter=ts_counter)
    elseif model_id === :level_1
        return DropletParams(T_inf=T_inf_C, U_rel=U_ms, Y_inf=y_inf, P_atm=P_ATM,
                             profile=:parabolico, as1989_blowing=false, model_variant=:level_1,
                             ts_solver=ts_solver, ts_counter=ts_counter)
    elseif model_id === :uniform
        return DropletParams(T_inf=T_inf_C, U_rel=U_ms, Y_inf=y_inf, P_atm=P_ATM,
                             profile=:uniforme, as1989_blowing=true, model_variant=:abramzon_sirignano,
                             ts_solver=:uniforme, ts_counter=ts_counter)
    elseif model_id === :iso_twb
        return DropletParams(T_inf=T_inf_C, U_rel=U_ms, Y_inf=y_inf, P_atm=P_ATM,
                             profile=:uniforme, as1989_blowing=true, model_variant=:abramzon_sirignano,
                             isothermal=true, T_d_fixed=twb, ts_solver=:uniforme, ts_counter=ts_counter)
    else
        error("Unknown physical model_id: $model_id")
    end
end

"""
    solve_jaiio_trajectory(model_id, d0_m, U_ms, T0_C, T_inf_C, HR; ...) -> NamedTuple

Executes a single ODE/DAE/Operator solve and returns the raw solution object without postprocessing,
ideal for accurate wall-time benchmarking.
"""
function solve_jaiio_trajectory(
    model_id::Symbol, d0_m::Float64, U_ms::Float64, T0_C::Float64, T_inf_C::Float64, HR::Float64;
    t_max::Float64=30000.0,
    ts_solver::Symbol=:brent,
    ts_counter::Union{Nothing, TsSolverCounter}=nothing,
    surrogate_coeffs=(modelo_6b_refit=nothing, modelo_7=nothing),
    pinn_bundle::Union{Nothing, PINNSurrogateBundle}=nothing,
    rtol::Float64=1e-6,
    atol::Float64=1e-9
)
    twb = stull_twb(T_inf_C, HR * 100.0)
    cat = model_category(model_id)

    if cat === :physical
        if ts_counter !== nothing
            reset!(ts_counter)
        end
        params = build_physical_params(model_id, d0_m, U_ms, T0_C, T_inf_C, HR;
                                       ts_solver=ts_solver, ts_counter=ts_counter)
        rho_l0 = liquid_density_water((model_id === :iso_twb ? twb : T0_C) + T_ZERO_C)
        m0     = rho_l0 * (π / 6.0) * (d0_m^3)
        T_init = model_id === :iso_twb ? twb : T0_C
        ode_alg = (ts_solver === :dae) ? nothing : AutoTsit5(Rosenbrock23(autodiff=false))
        ode_atol = (ts_solver === :dae) ? atol : [min(1e-18, m0 * 1e-9), atol]
        sol    = solve_droplet(params, [m0, T_init], (0.0, t_max);
                               solver=ode_alg, rtol=rtol, atol=ode_atol, maxiters=10^6)
        return (sol=sol, params=params, m0=m0, twb=twb, category=:physical)
    elseif cat === :neural
        bundle = pinn_bundle === nothing ? load_pinn_models() : pinn_bundle
        d0_mm = d0_m * 1e3
        if model_id === :pinn_c
            sol = solve_pinn_c(bundle.pinn_c, d0_mm, U_ms, T0_C, T_inf_C, HR;
                               t_end=t_max, solver=Tsit5(), rtol=rtol, atol=atol)
            return (sol=sol, params=bundle.pinn_c, m0=(π / 6.0) * (d0_mm^3), twb=twb, category=:neural)
        elseif model_id === :pinn_2eq
            sol = solve_pinn_2eq(bundle.pinn_2eq, d0_mm, U_ms, T0_C, T_inf_C, HR;
                                 t_end=t_max, solver=AutoTsit5(Rosenbrock23(autodiff=false)), rtol=rtol, atol=atol)
            return (sol=sol, params=bundle.pinn_2eq, m0=(π / 6.0) * (d0_mm^3), twb=twb, category=:neural)
        elseif model_id === :pi_deeponet
            op_res = deeponet_eval_case(bundle.pi_deeponet, d0_mm, U_ms, T0_C, T_inf_C, HR)
            sol = (t=[0.0, op_res.t_evap_s], nfev=1, op=op_res)
            return (sol=sol, params=bundle.pi_deeponet, m0=(π / 6.0) * (d0_mm^3), twb=twb, category=:neural)
        elseif model_id === :resmlp
            op_res = resmlp_eval_case(bundle.resmlp, d0_mm, U_ms, T0_C, T_inf_C, HR)
            sol = (t=[0.0, op_res.t_evap_s], nfev=1, op=op_res)
            return (sol=sol, params=bundle.resmlp, m0=(π / 6.0) * (d0_mm^3), twb=twb, category=:neural)
        else
            error("Unknown neural model_id: $model_id")
        end
    else
        c_override = model_id === :modelo_6b_refit ? surrogate_coeffs.modelo_6b_refit :
                     model_id === :modelo_7        ? surrogate_coeffs.modelo_7        : nothing
        cfg = SurrogateConfig(
            model_type = model_id,
            d0_mm      = d0_m * 1e3,
            T_inf      = T_inf_C,
            T_wb       = twb,
            U          = U_ms,
            HR         = HR,
            T0         = T0_C,
            t_end      = t_max,
            coeffs     = c_override,
        )
        sol = solve_surrogate(cfg; solver=Tsit5(), rtol=rtol, atol=atol)
        return (sol=sol, params=cfg, m0=(π / 6.0) * ((d0_m * 1e3)^3), twb=twb, category=:surrogate)
    end
end

"""
    run_jaiio_case(model_id, d0_m, U_ms, T0_C, T_inf_C, HR; ...) -> NamedTuple

Runs a case, times the raw ODE/Operator solve, and computes all t=0 and integrated summary metrics.
"""
function run_jaiio_case(
    model_id::Symbol, d0_m::Float64, U_ms::Float64, T0_C::Float64, T_inf_C::Float64, HR::Float64;
    t_max::Float64=30000.0,
    ts_solver::Symbol=:brent,
    count_ts::Bool=false,
    surrogate_coeffs=(modelo_6b_refit=nothing, modelo_7=nothing),
    pinn_bundle::Union{Nothing, PINNSurrogateBundle}=nothing,
    rtol::Float64=1e-6,
    atol::Float64=1e-9
)
    ctr = count_ts ? TsSolverCounter() : nothing
    t0_ns = time_ns()
    res = solve_jaiio_trajectory(
        model_id, d0_m, U_ms, T0_C, T_inf_C, HR;
        t_max=t_max, ts_solver=ts_solver, ts_counter=ctr,
        surrogate_coeffs=surrogate_coeffs, pinn_bundle=pinn_bundle, rtol=rtol, atol=atol
    )
    elapsed_s = (time_ns() - t0_ns) / 1e9

    sol = res.sol
    t_evap = Float64(sol.t[end])
    nfev   = _extract_nfev(sol)
    n_steps = length(sol.t)
    n_ts_calls = ctr === nothing ? 0 : ctr.calls
    n_ts_iters = ctr === nothing ? 0 : ctr.iters
    total_ops  = nfev + n_ts_iters
    cwu        = total_ops * cwu_weight(model_id)
    ctr !== nothing && reset!(ctr)

    V0_mm3 = (π / 6.0) * ((d0_m * 1e3)^3)

    if res.category === :physical
        p = res.params
        m0 = res.m0
        T_init = model_id === :iso_twb ? res.twb : T0_C
        Ts0 = if model_id === :iso_twb
            res.twb
        elseif model_id === :uniform
            T0_C
        elseif p.ts_solver === :dae && length(sol.u[1]) >= 3
            sol.u[1][3]
        else
            Float64(solve_surface_temperature(m0, T_init, p))
        end
        alg0 = calc_droplet_algebraic(m0, T_init, Ts0, p)

        m_end = max(Float64(sol.u[end][1]), 0.0)
        T_end = Float64(sol.u[end][2])
        rl_end = liquid_density_water(T_end + T_ZERO_C)
        V_end_mm3 = (m_end / rl_end) * 1e9
        dV_dt_bulk = t_evap > 0.0 ? (V0_mm3 - V_end_mm3) / t_evap : 0.0

        mdots = Float64[]
        sizehint!(mdots, n_steps)
        for i in 1:n_steps
            mi = sol.u[i][1]
            Ti = sol.u[i][2]
            if mi > 1e-18
                Tsi = if model_id === :iso_twb
                    res.twb
                elseif model_id === :uniform
                    Ti
                elseif p.ts_solver === :dae && length(sol.u[i]) >= 3
                    sol.u[i][3]
                else
                    Float64(solve_surface_temperature(mi, Ti, p))
                end
                push!(mdots, Float64(calc_droplet_algebraic(mi, Ti, Tsi, p).m_dot))
            end
        end
        mdot_med = isempty(mdots) ? alg0.m_dot : median(mdots)

        return (
            label = string(model_id),
            d0_m = d0_m, U_ms = U_ms, T0_C = T0_C, T_inf_C = T_inf_C, HR = HR, T_wb_C = res.twb,
            t_evap_s = t_evap,
            completed = (t_evap < t_max - 1.0),
            m_dot_0_kgs = Float64(alg0.m_dot),
            m_dot_median_kgs = Float64(mdot_med),
            dV_dt_bulk_mm3s = Float64(dV_dt_bulk),
            T_s0_C = Float64(Ts0),
            B_M_0 = Float64(alg0.B_M),
            B_T_0 = Float64(alg0.B_T),
            Re_0 = Float64(alg0.Re),
            Nu_star_0 = Float64(alg0.Nu_star),
            Sh_star_0 = Float64(alg0.Sh_star),
            q_conv_0_W = Float64(alg0.q_g),
            q_lat_0_W = Float64(alg0.m_dot * alg0.DHev),
            elapsed_s = elapsed_s,
            nfev = nfev,
            n_timesteps = n_steps,
            n_ts_calls = n_ts_calls,
            n_ts_iters = n_ts_iters,
            total_ops = total_ops,
            cwu = cwu,
        )
    elseif res.category === :neural
        Re0 = Float64(calc_reynolds_surrogate(d0_m * 1e3, U_ms, T_inf_C))
        d0_mm = d0_m * 1e3
        sc = PINNSurrogates.compute_static_thermo_scales(T_inf_C, HR, T0_C)

        if model_id === :pinn_c
            m_c = res.params
            rho_l0 = liquid_density_water(res.twb + T_ZERO_C)
            C0 = PINNSurrogates.pinn_c_eval_C(m_c, d0_mm, d0_mm, U_ms, T0_C, T_inf_C, HR, sc.T_wb, sc.Theta0, sc.BM_wb)
            dV0 = -C0 * sqrt(sphere_surface_from_d(d0_mm)) * (T_inf_C - sc.T_wb)
            m_dot_0 = -dV0 * 1e-9 * rho_l0
            mdots = Float64[]
            for i in 1:n_steps
                Vi = sol.u[i][1]
                if Vi > 1e-9
                    di = sphere_diameter_from_volume(Vi)
                    Ci = PINNSurrogates.pinn_c_eval_C(m_c, di, d0_mm, U_ms, T0_C, T_inf_C, HR, sc.T_wb, sc.Theta0, sc.BM_wb)
                    dVi = -Ci * sqrt(sphere_surface_from_d(di)) * (T_inf_C - sc.T_wb)
                    push!(mdots, -dVi * 1e-9 * rho_l0)
                end
            end
            mdot_med = isempty(mdots) ? m_dot_0 : median(mdots)
            V_end_mm3 = max(Float64(sol.u[end][1]), 0.0)
            dV_dt_bulk = t_evap > 0.0 ? (V0_mm3 - V_end_mm3) / t_evap : 0.0
            Ts0 = res.twb
        elseif model_id === :pinn_2eq
            m_2 = res.params
            st0 = pinn_2eq_step_physics(m_2, V0_mm3, T0_C, d0_mm, U_ms, T_inf_C, HR, sc.T_wb, sc.Y_inf)
            Ts0 = Float64(st0.Ts_hat)
            m_dot_0 = Float64(st0.m_dot)
            mdots = Float64[]
            for i in 1:n_steps
                Vi, Ti = sol.u[i][1], sol.u[i][2]
                if Vi > 1e-9
                    sti = pinn_2eq_step_physics(m_2, Vi, Ti, d0_mm, U_ms, T_inf_C, HR, sc.T_wb, sc.Y_inf)
                    push!(mdots, Float64(sti.m_dot))
                end
            end
            mdot_med = isempty(mdots) ? m_dot_0 : median(mdots)
            V_end_mm3 = max(Float64(sol.u[end][1]), 0.0)
            dV_dt_bulk = t_evap > 0.0 ? (V0_mm3 - V_end_mm3) / t_evap : 0.0
        elseif model_id === :pi_deeponet
            op = sol.op
            Ts0 = Float64(op.T_s0_C)
            m_dot_0 = Float64(op.m_dot_0_kgs)
            dV_dt_bulk = Float64(op.dV_dt_bulk_mm3s)
            mdot_med = dV_dt_bulk * 1e-9 * liquid_density_water(res.twb + T_ZERO_C)
        else
            op = sol.op
            Ts0 = Float64(op.T_s0_C)
            m_dot_0 = Float64(op.m_dot_0_kgs)
            dV_dt_bulk = Float64(op.dV_dt_bulk_mm3s)
            mdot_med = dV_dt_bulk * 1e-9 * liquid_density_water(res.twb + T_ZERO_C)
        end

        return (
            label = string(model_id),
            d0_m = d0_m, U_ms = U_ms, T0_C = T0_C, T_inf_C = T_inf_C, HR = HR, T_wb_C = res.twb,
            t_evap_s = t_evap,
            completed = (t_evap < t_max - 1.0),
            m_dot_0_kgs = Float64(m_dot_0),
            m_dot_median_kgs = Float64(mdot_med),
            dV_dt_bulk_mm3s = Float64(dV_dt_bulk),
            T_s0_C = Float64(Ts0),
            B_M_0 = NaN,
            B_T_0 = NaN,
            Re_0 = Re0,
            Nu_star_0 = NaN,
            Sh_star_0 = NaN,
            q_conv_0_W = NaN,
            q_lat_0_W = NaN,
            elapsed_s = elapsed_s,
            nfev = nfev,
            n_timesteps = n_steps,
            n_ts_calls = 0,
            n_ts_iters = 0,
            total_ops = total_ops,
            cwu = cwu,
        )
    else
        cfg = res.params
        du = zeros(length(sol.u[1]))
        Surrogates.surrogate_ode_rhs!(du, sol.u[1], cfg, 0.0)
        rho_l0 = liquid_density_water(res.twb + T_ZERO_C)
        m_dot_0 = -du[1] * 1e-9 * rho_l0

        mdots = Float64[]
        sizehint!(mdots, n_steps)
        for i in 1:n_steps
            if sol.u[i][1] > 1e-9
                Surrogates.surrogate_ode_rhs!(du, sol.u[i], cfg, sol.t[i])
                push!(mdots, -du[1] * 1e-9 * rho_l0)
            end
        end
        mdot_med = isempty(mdots) ? m_dot_0 : median(mdots)
        V_end_mm3 = max(Float64(sol.u[end][1]), 0.0)
        dV_dt_bulk = t_evap > 0.0 ? (V0_mm3 - V_end_mm3) / t_evap : 0.0
        Re0 = Float64(calc_reynolds_surrogate(d0_m * 1e3, U_ms, T_inf_C))
        Ts0 = res.twb

        return (
            label = string(model_id),
            d0_m = d0_m, U_ms = U_ms, T0_C = T0_C, T_inf_C = T_inf_C, HR = HR, T_wb_C = res.twb,
            t_evap_s = t_evap,
            completed = (t_evap < t_max - 1.0),
            m_dot_0_kgs = Float64(m_dot_0),
            m_dot_median_kgs = Float64(mdot_med),
            dV_dt_bulk_mm3s = Float64(dV_dt_bulk),
            T_s0_C = Float64(Ts0),
            B_M_0 = NaN,
            B_T_0 = NaN,
            Re_0 = Re0,
            Nu_star_0 = NaN,
            Sh_star_0 = NaN,
            q_conv_0_W = NaN,
            q_lat_0_W = NaN,
            elapsed_s = elapsed_s,
            nfev = nfev,
            n_timesteps = n_steps,
            n_ts_calls = 0,
            n_ts_iters = 0,
            total_ops = total_ops,
            cwu = cwu,
        )
    end
end

end # module JAIIOCommon
