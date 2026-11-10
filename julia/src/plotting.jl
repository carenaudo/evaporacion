"""
    Evaporacion.Plotting

Plotting utilities using Plots.jl for evaporation trajectories, profiles, and parity charts.
"""
module Plotting

using Plots

export plot_droplet_trajectory, plot_sessile_evolution, plot_temperature_profile

"""
    plot_droplet_trajectory(df; savepath=nothing)

Plots d^2(t), mean droplet temperature T_bar(t), and surface temperature Ts(t).
"""
function plot_droplet_trajectory(df; savepath=nothing)
    d2_mm2 = (df.d .* 1000.0) .^ 2
    p1 = plot(df.t, d2_mm2, xlabel="Time [s]", ylabel="d² [mm²]", label="d²(t)", lw=2, color=:blue)
    p2 = plot(df.t, df.T_bar, xlabel="Time [s]", ylabel="Temperature [°C]", label="T_bar", lw=2, color=:red)
    plot!(p2, df.t, df.Ts, label="Ts", lw=2, ls=:dash, color=:orange)

    p = plot(p1, p2, layout=(2, 1), size=(800, 600))
    if savepath !== nothing
        savefig(p, savepath)
    end
    return p
end

"""
    plot_sessile_evolution(sol; savepath=nothing)

Plots contact radius R(t) and contact angle θ(t) for sessile droplets.
"""
function plot_sessile_evolution(sol; savepath=nothing)
    p1 = plot(sol.t, sol.R .* 1000.0, xlabel="Time [s]", ylabel="Contact Radius R [mm]", lw=2, color=:blue, label="R(t)")
    p2 = plot(sol.t, rad2deg.(sol.theta), xlabel="Time [s]", ylabel="Contact Angle θ [deg]", lw=2, color=:green, label="θ(t)")

    p = plot(p1, p2, layout=(2, 1), size=(800, 600), title="Sessile Evaporation ($(sol.mode) Mode)")
    if savepath !== nothing
        savefig(p, savepath)
    end
    return p
end

"""
    plot_temperature_profile(sol_pde; savepath=nothing)

Plots internal temperature profile T(ξ, t) across radial coordinates.
"""
function plot_temperature_profile(sol_pde; savepath=nothing)
    n_times, n_nodes = size(sol_pde.T_profile)
    xi = range(0.0, 1.0, length=n_nodes)
    p = plot(xlabel="Normalized radius ξ = r / R(t)", ylabel="Temperature [°C]", title="Internal Temperature Profiles", legend=:topleft)

    step = max(1, n_times ÷ 6)
    for i in 1:step:n_times
        plot!(p, xi, sol_pde.T_profile[i, :], label="t = $(round(sol_pde.t[i], digits=1)) s", lw=1.5)
    end

    if savepath !== nothing
        savefig(p, savepath)
    end
    return p
end

end # module Plotting
