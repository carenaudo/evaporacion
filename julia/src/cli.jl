"""
    Evaporacion.CLI

Command-line interface for the Evaporacion.jl suite.
"""
module CLI

using ArgParse
using CSV
using DataFrames
using ..Common
using ..Droplet
using ..Surrogates
using ..Sessile
using ..RadialPDE
using ..Plotting

export parse_commandline, main

function parse_commandline(args=ARGS)
    s = ArgParseSettings(
        description = "Evaporacion.jl — Advanced Droplet Evaporation Simulation Suite in High-Performance Julia"
    )

    @add_arg_table! s begin
        "--model", "-m"
            help = "Model type: droplet, surrogate, sessile, radial_pde"
            default = "droplet"
        "--d0"
            help = "Initial droplet diameter in mm"
            arg_type = Float64
            default = 1.0
        "--Tinf"
            help = "Ambient temperature in °C"
            arg_type = Float64
            default = 35.0
        "--T0"
            help = "Initial droplet temperature in °C"
            arg_type = Float64
            default = 22.0
        "--U"
            help = "Relative velocity in m/s"
            arg_type = Float64
            default = 0.5
        "--HR"
            help = "Relative humidity (percentage 0-100 or fraction 0-1)"
            arg_type = Float64
            default = 40.0
        "--profile"
            help = "Droplet internal profile: uniforme, parabolico, cubico"
            default = "parabolico"
        "--surrogate-type"
            help = "Surrogate model: modelo_2b, modelo_4, modelo_5, modelo_6, modelo_6b_refit, modelo_7"
            default = "modelo_6b_refit"
        "--sessile-mode"
            help = "Sessile evaporation mode: CR, CA, SS, SJ"
            default = "CR"
        "--t-end"
            help = "Maximum simulation time in seconds"
            arg_type = Float64
            default = 3000.0
        "--save-csv"
            help = "Path to export output trajectory as CSV"
            arg_type = String
            default = ""
        "--plot"
            help = "Path to export PNG visualization"
            arg_type = String
            default = ""
    end

    return parse_args(args, s)
end

function main(args=ARGS)
    parsed = parse_commandline(args)

    hr_frac = parsed["HR"] > 1.0 ? parsed["HR"] / 100.0 : parsed["HR"]
    d0_m = parsed["d0"] / 1000.0

    if parsed["model"] == "droplet"
        println("[INFO] Running Rigorous Spherical Droplet Model...")
        p_sat_inf = psat_water(parsed["Tinf"] + 273.15)
        Y_inf = surface_mass_fraction_from_psat(hr_frac * p_sat_inf)

        prof_sym = Symbol(parsed["profile"])
        params = DropletParams(
            T_inf = parsed["Tinf"],
            U_rel = parsed["U"],
            Y_inf = Y_inf,
            profile = prof_sym
        )

        rho_l0 = liquid_density_water(parsed["T0"] + 273.15)
        m0 = rho_l0 * (π / 6.0) * (d0_m^3)
        u0 = [m0, parsed["T0"]]

        sol = solve_droplet(params, u0, (0.0, parsed["t-end"]))
        df = postprocess_trajectory(sol, params)
        println("[SUCCESS] Integration complete: $(nrow(df)) steps, final time = $(round(df.t[end], digits=2)) s")

        if !isempty(parsed["save-csv"])
            CSV.write(parsed["save-csv"], df)
            println("[INFO] Saved results to $(parsed["save-csv"])")
        end
        if !isempty(parsed["plot"])
            plot_droplet_trajectory(df; savepath=parsed["plot"])
            println("[INFO] Saved plot to $(parsed["plot"])")
        end

    elseif parsed["model"] == "surrogate"
        println("[INFO] Running Surrogate Model ($(parsed["surrogate-type"]))...")
        T_wb = stull_twb(parsed["Tinf"], hr_frac * 100.0)
        cfg = SurrogateConfig(
            model_type = Symbol(parsed["surrogate-type"]),
            d0_mm = parsed["d0"],
            T_inf = parsed["Tinf"],
            T_wb = T_wb,
            U = parsed["U"],
            HR = hr_frac,
            T0 = parsed["T0"],
            t_end = parsed["t-end"]
        )
        sol = solve_surrogate(cfg)
        println("[SUCCESS] Surrogate integration complete: final time = $(round(sol.t[end], digits=2)) s")

    elseif parsed["model"] == "sessile"
        println("[INFO] Running Sessile Droplet Model ($(parsed["sessile-mode"]))...")
        cfg = SessileDropletConfig(
            R0 = d0_m / 2.0,
            theta0 = 0.5 * π,
            mode = Symbol(parsed["sessile-mode"]),
            t_end = parsed["t-end"]
        )
        sol = solve_sessile(cfg)
        println("[SUCCESS] Sessile integration complete: final time = $(round(sol.t[end], digits=2)) s, analytical t_evap = $(round(sol.t_evap_analytical, digits=2)) s")

    elseif parsed["model"] == "radial_pde"
        println("[INFO] Running 1D Moving-Boundary Radial PDE Model...")
        cfg = RadialPDEConfig(
            d0 = d0_m,
            T_drop_0 = parsed["T0"],
            T_air = parsed["Tinf"],
            RH = hr_frac,
            t_end = parsed["t-end"]
        )
        sol = solve_radial_pde(cfg)
        println("[SUCCESS] Radial PDE integration complete: final time = $(round(sol.t[end], digits=2)) s")
    end
end

end # module CLI
