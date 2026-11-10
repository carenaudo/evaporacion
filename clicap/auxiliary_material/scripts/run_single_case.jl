"""
Single Case Runner - Auxiliary Material Script
Allows running any single CLICAP environmental case (or custom parameters) with any implemented model in Evaporacion.jl.

Usage examples:
  julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain inflight --case 1 --d0 200 --model abramzon_sirignano
  julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --domain sessile --case 7 --theta0 40 --model CR
  julia --project=julia clicap/auxiliary_material/scripts/run_single_case.jl --T 25.0 --HR 50.0 --d0 600 --model d2_convective
"""

using ArgParse
using Printf
using DataFrames
using CSV

# Load Evaporacion.jl package
include(joinpath(@__DIR__, "..", "..", "..", "julia", "src", "Evaporacion.jl"))
using .Evaporacion

# Standard 9 Environmental Cases from Renaudo et al. (Table 1)
const CLICAP_CASES = Dict(
    1 => (T_K = 283.15, HR = 20.0),
    2 => (T_K = 283.15, HR = 50.0),
    3 => (T_K = 283.15, HR = 80.0),
    4 => (T_K = 298.15, HR = 20.0),
    5 => (T_K = 298.15, HR = 50.0),
    6 => (T_K = 298.15, HR = 80.0),
    7 => (T_K = 313.15, HR = 20.0),
    8 => (T_K = 313.15, HR = 50.0),
    9 => (T_K = 313.15, HR = 80.0)
)

function parse_cli_args()
    s = ArgParseSettings(
        description = "Run a single evaporation simulation case across in-flight or sessile droplet models."
    )

    @add_arg_table! s begin
        "--domain"
            help = "Domain of evaporation: 'inflight' or 'sessile'"
            arg_type = String
            default = "inflight"
        "--case"
            help = "CLICAP case index (1 to 9). Overrides --T and --HR if provided."
            arg_type = Int
            default = 0
        "--T"
            help = "Ambient air temperature in Celsius (used if --case 0)"
            arg_type = Float64
            default = 25.0
        "--HR"
            help = "Ambient relative humidity in percent (0 to 100, used if --case 0)"
            arg_type = Float64
            default = 50.0
        "--P_atm"
            help = "Atmospheric pressure in Pa"
            arg_type = Float64
            default = 101325.0
        "--d0"
            help = "Initial droplet diameter in micrometers (for inflight or equivalent spherical volume for sessile)"
            arg_type = Float64
            default = 200.0
        "--theta0"
            help = "Initial contact angle in degrees for sessile droplets (e.g. 20, 40, 80, 90)"
            arg_type = Float64
            default = 40.0
        "--model"
            help = "Model identifier:\n" *
                   "  In-flight: d2_classic, d2_stefan, d2_convective, d2_revised, abramzon_sirignano, abramzon_parabolic, zhou_cubic, vlasov_kinetic, surrogate_6b, surrogate_8dyn\n" *
                   "  Sessile: CA, CR, SS, SJ, thin_cr, thin_ca, picknett, substrate_cooling, vapor_shielding, saxton_tanner"
            arg_type = String
            default = "abramzon_sirignano"
        "--U_rel"
            help = "Relative air velocity in m/s (if 0 and model has convection, terminal fall velocity is computed dynamically)"
            arg_type = Float64
            default = 0.0
        "--output"
            help = "Optional CSV file path to export result"
            arg_type = String
            default = ""
    end

    return parse_args(s)
end

function run_inflight(T_C, HR_pct, P_atm, d0_um, model_name, U_rel)
    d0_m = d0_um * 1e-6
    T_K = T_C + 273.15
    HR_frac = HR_pct / 100.0
    psat_inf = psat_water(T_K)
    pv_inf = HR_frac * psat_inf
    M_w = Evaporacion.Common.Constants.M_WATER
    M_a = Evaporacion.Common.Constants.M_AIR
    Y_inf = (pv_inf * M_w) / (pv_inf * M_w + (P_atm - pv_inf) * M_a)
    T_wb_C = stull_twb(T_C, HR_pct)

    println("================================================================================")
    println("RUNNING IN-FLIGHT DROPLET CASE")
    println("================================================================================")
    @printf("Ambient Conditions : T = %.2f °C (%.2f K), HR = %.1f%%, P = %.0f Pa\n", T_C, T_K, HR_pct, P_atm)
    @printf("Psychrometrics     : T_wb = %.2f °C, Y_inf = %.5f kg/kg\n", T_wb_C, Y_inf)
    @printf("Initial Droplet    : d0 = %.1f µm (%.3e m)\n", d0_um, d0_m)
    @printf("Model Selection    : %s\n", model_name)
    println("--------------------------------------------------------------------------------")

    t_evap = 0.0
    K_initial = 0.0

    if model_name in ["d2_classic", "d2_stefan", "d2_convective", "d2_revised", "abramzon_sirignano", "abramzon_parabolic", "zhou_cubic"]
        m_sym = if model_name == "d2_classic"
            :d2_classic
        elseif model_name == "d2_stefan"
            :d2_stefan
        elseif model_name == "d2_convective"
            :d2_convective
        elseif model_name == "d2_revised"
            :d2_revised_dallabarba
        else
            :abramzon_sirignano
        end

        profile_sym = if model_name == "abramzon_uniform"
            :uniforme
        elseif model_name == "zhou_cubic"
            :zhou_cubic
        else
            :parabolico
        end

        isothermal = m_sym in (:d2_classic, :d2_stefan, :d2_convective, :d2_revised_dallabarba)
        use_terminal = (U_rel == 0.0) && (m_sym in (:d2_convective, :d2_revised_dallabarba, :abramzon_sirignano))
        T_init_C = isothermal ? T_wb_C : T_C

        rho_l0 = liquid_density_water(T_init_C + 273.15)
        m0 = rho_l0 * (π / 6.0) * (d0_m^3)

        params = DropletParams(
            T_inf = T_C,
            Y_inf = Y_inf,
            P_atm = P_atm,
            U_rel = U_rel,
            model_variant = m_sym,
            profile = profile_sym,
            terminal_velocity = use_terminal,
            isothermal = isothermal,
            T_d_fixed = T_wb_C
        )

        # Initial algebraic rate
        alg0 = calc_droplet_algebraic(m0, T_init_C, T_init_C, params)
        K_initial = (4.0 * alg0.m_dot / (π * alg0.rho_l * alg0.d)) * 1e6 # mm^2/s

        t_est = (d0_um / 50.0)^2 * 10.0
        tspan = (0.0, max(50.0, t_est * 5.0))

        sol = solve_droplet(params, [m0, T_init_C], tspan)
        t_evap = sol.t[end]

    elseif model_name == "vlasov_kinetic"
        # Mixed kinetic-diffusive model (Vlasov 2021)
        rho_l0 = liquid_density_water(T_wb_C + 273.15)
        D_val = gas_diffusivity_water_air(T_K, P_atm)
        delta_c = (psat_water(T_wb_C + 273.15) / (Evaporacion.Common.Constants.R_GAS * (T_wb_C + 273.15)) * M_w) * (1.0 - HR_frac)
        R0 = d0_m / 2.0
        alpha_evap = 0.04
        v_th = sqrt(8.0 * Evaporacion.Common.Constants.R_GAS * T_K / (π * M_w))

        # Analytical lifetime extension from Vlasov
        tau_diff = (rho_l0 * R0^2) / (2.0 * D_val * delta_c)
        l_k = vlasov_kinetic_length(D_val, alpha_evap, v_th)
        t_evap = tau_diff * (1.0 + 2.0 * l_k / R0)
        K_initial = (d0_m * 1e3)^2 / t_evap

    elseif model_name in ["surrogate_6b", "surrogate_8dyn"]
        stype = model_name == "surrogate_6b" ? :modelo_6b_refit : :modelo_8_dynamic
        cfg = SurrogateConfig(
            model_type = stype,
            d0_mm = d0_um * 1e-3,
            T_inf = T_C,
            T_wb = T_wb_C,
            U = max(0.2, U_rel),
            HR = HR_frac,
            T0 = T_C,
            t_end = 2000.0
        )
        sol = solve_surrogate(cfg)
        t_evap = sol.t[end]
        K_initial = (d0_um * 1e-3)^2 / t_evap
    else
        error("Unknown in-flight model: $model_name")
    end

    # Also compute in-flight extinction parameters (Renaudo 2022, 2023)
    k_ext = renaudo_evaporation_parameter(T_C, HR_pct)
    d_min_50cm = renaudo_extinction_diameter(0.5, T_C, HR_pct) * 1e6 # µm
    d_min_100cm = renaudo_extinction_diameter(1.0, T_C, HR_pct) * 1e6 # µm
    t_dep = renaudo_deposition_time(d0_m, 0.5, T_C, HR_pct)
    d_land_um = renaudo_landing_diameter(d0_m, 0.5, T_C, HR_pct) * 1e6

    @printf("SIMULATION RESULTS:\n")
    @printf("  Evaporation Lifetime (t_evap)       : %.3f s\n", t_evap)
    @printf("  Initial Evaporation Rate K(t=0)     : %.5e mm²/s\n", K_initial)
    @printf("  Average Evaporation Rate K_avg      : %.5e mm²/s\n", (d0_um*1e-3)^2 / t_evap)
    println("--------------------------------------------------------------------------------")
    @printf("ATMOSPHERIC EXTINCTION (Renaudo 2022):\n")
    @printf("  Evaporation Parameter k             : %.3e s/m²\n", k_ext)
    @printf("  Cut-off Diameter d_min (H = 0.5 m)  : %.2f µm\n", d_min_50cm)
    @printf("  Cut-off Diameter d_min (H = 1.0 m)  : %.2f µm\n", d_min_100cm)
    @printf("  Sedimentation Landing Time (H=0.5m) : %s\n", isinf(t_dep) ? "Never (evaporates in flight)" : @sprintf("%.3f s", t_dep))
    @printf("  Final Soil Landing Diameter         : %.2f µm (%.1f%% of initial)\n", d_land_um, (d_land_um/d0_um)*100.0)
    println("================================================================================")

    return (
        domain = "inflight",
        model = model_name,
        T_C = T_C,
        HR_pct = HR_pct,
        d0_um = d0_um,
        t_evap_s = t_evap,
        K_t0_mm2_s = K_initial,
        d_min_50cm_um = d_min_50cm,
        d_land_um = d_land_um
    )
end

function run_sessile(T_C, HR_pct, P_atm, d0_equiv_um, theta0_deg, model_name)
    d0_m = d0_equiv_um * 1e-6
    V0 = (π / 6.0) * (d0_m^3)
    T_K = T_C + 273.15
    HR_frac = HR_pct / 100.0

    rho_l = liquid_density_water(T_K)
    D_val = gas_diffusivity_water_air(T_K, P_atm)
    psat = psat_water(T_K)
    c_sat = (psat * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * T_K)
    delta_c = (1.0 - HR_frac) * c_sat

    th0_rad = deg2rad(theta0_deg)
    R0 = r_from_volume_sessile(V0, th0_rad)

    println("================================================================================")
    println("RUNNING SESSILE DROPLET CASE")
    println("================================================================================")
    @printf("Ambient Conditions : T = %.2f °C, HR = %.1f%%, P = %.0f Pa\n", T_C, HR_pct, P_atm)
    @printf("Droplet Geometry   : V0 = %.3e m³ (equiv d0 = %.1f µm)\n", V0, d0_equiv_um)
    @printf("Contact Line       : θ₀ = %.1f° (%.4f rad), R₀ = %.2f µm\n", theta0_deg, th0_rad, R0*1e6)
    @printf("Model Selection    : %s\n", model_name)
    println("--------------------------------------------------------------------------------")

    t_evap = 0.0

    if model_name == "CA"
        cfg = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CA, t_end=2000.0)
        sol = solve_sessile(cfg)
        t_evap = sol.t_evap_analytical
    elseif model_name == "CR"
        cfg = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CR, t_end=2000.0)
        sol = solve_sessile(cfg)
        t_evap = sol.t_evap_analytical
    elseif model_name == "SS"
        cfg = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, theta_rec=th0_rad/2.0, mode=:SS, t_end=2000.0)
        sol = solve_sessile(cfg)
        t_evap = sol.t_evap_analytical
    elseif model_name == "SJ"
        cfg = SessileDropletConfig(D=D_val, delta_c=delta_c, rho_l=rho_l, R0=R0, theta0=th0_rad, theta_min=th0_rad/2.0, theta_max=th0_rad, mode=:SJ, t_end=2000.0)
        sol = solve_sessile(cfg)
        t_evap = sol.t_evap_analytical
    elseif model_name == "thin_cr"
        t_evap = lifetime_cr_thin(R0, th0_rad, D_val, delta_c, rho_l)
    elseif model_name == "thin_ca"
        t_evap = lifetime_ca_thin(R0, th0_rad, D_val, delta_c, rho_l)
    elseif model_name == "picknett"
        # Capacitance polynomial model
        C_ratio = picknett_bexon_capacitance(th0_rad)
        t_cr_est = (rho_l * V0) / (4.0 * π * R0 * D_val * delta_c * C_ratio)
        t_evap = t_cr_est
    elseif model_name == "substrate_cooling"
        # Substrate cooling on glass (h_s = 1 mm, k_s = 1.4 W/mK)
        L_val = latent_heat_water(T_K)
        csat_plus = (psat_water(T_K + 0.5) * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * (T_K + 0.5))
        csat_minus = (psat_water(T_K - 0.5) * Evaporacion.Common.Constants.M_WATER) / (Evaporacion.Common.Constants.R_GAS * (T_K - 0.5))
        dcsat_dT_val = csat_plus - csat_minus
        tparams = SessileThermalParams(L=L_val, h_s=1e-3, k_s_eff=1.4, dcsat_dT=dcsat_dT_val)
        delta_c_eff = cooling_effective_delta_c(R0, th0_rad, D_val, delta_c, tparams)
        cfg_cool = SessileDropletConfig(D=D_val, delta_c=delta_c_eff, rho_l=rho_l, R0=R0, theta0=th0_rad, mode=:CR, t_end=2000.0)
        sol_cool = solve_sessile(cfg_cool)
        t_evap = sol_cool.t_evap_analytical
    elseif model_name == "vapor_shielding"
        # Shielding for a pair of droplets at separation B = 3*R0
        t_evap = lifetime_pair_cr(3.0 * R0, R0, th0_rad, D_val, delta_c, rho_l)
    else
        error("Unknown sessile model: $model_name")
    end

    @printf("SIMULATION RESULTS:\n")
    @printf("  Evaporation Lifetime (t_evap)       : %.3f s\n", t_evap)
    println("================================================================================")

    return (
        domain = "sessile",
        model = model_name,
        T_C = T_C,
        HR_pct = HR_pct,
        d0_equiv_um = d0_equiv_um,
        theta0_deg = theta0_deg,
        t_evap_s = t_evap
    )
end

function main()
    args = parse_cli_args()

    # Determine ambient conditions
    case_idx = args["case"]
    T_C = args["T"]
    HR_pct = args["HR"]

    if case_idx in 1:9
        cdata = CLICAP_CASES[case_idx]
        T_C = cdata.T_K - 273.15
        HR_pct = cdata.HR
        println("Using predefined CLICAP Case $case_idx: T = $T_C °C ($(cdata.T_K) K), HR = $HR_pct%")
    end

    domain = lowercase(args["domain"])
    model_name = args["model"]
    d0_um = args["d0"]
    P_atm = args["P_atm"]

    res = if domain == "inflight"
        run_inflight(T_C, HR_pct, P_atm, d0_um, model_name, args["U_rel"])
    elseif domain == "sessile"
        run_sessile(T_C, HR_pct, P_atm, d0_um, args["theta0"], model_name)
    else
        error("Domain must be 'inflight' or 'sessile', got '$domain'")
    end

    out_path = args["output"]
    if !isempty(out_path)
        df = DataFrame([res])
        CSV.write(out_path, df)
        println("Exported results to: $out_path")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
