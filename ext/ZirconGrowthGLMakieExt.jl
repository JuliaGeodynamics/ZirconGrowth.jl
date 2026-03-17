module ZirconGrowthGLMakieExt

using ZirconGrowth
using GLMakie

# ─── Colour palette ───────────────────────────────────────────────────────────
_turbo_colors(n) = [Makie.to_color(GLMakie.Makie.ColorSchemes.turbo[i / n]) for i in 1:n]

# ─── Figure 1 equivalent: Kd and D vs 1000/T ─────────────────────────────────

function ZirconGrowth.plot_trace_properties(Th::AbstractVector{Float64},
                                            XH2O::Float64,
                                            element_names::NTuple{N,String};
                                            fig::Figure = Figure(size=(1200, 500))) where {N}
    nt = length(Th)
    xvals = 1000.0 ./ Th
    nelem = N

    kd_mat  = Matrix{Float64}(undef, nt, nelem)
    dif_mat = Matrix{Float64}(undef, nt, nelem)
    zr_kd   = Vector{Float64}(undef, nt)
    zr_dif  = Vector{Float64}(undef, nt)

    for i in 1:nt
        kd  = ZirconGrowth.partition_coefficients(Th[i])
        dif = ZirconGrowth.trace_diffusivities(Th[i], XH2O)
        Csat = ZirconGrowth.zr_saturation(Th[i])
        zr_kd[i]  = log(490_000.0 / Csat)
        zr_dif[i] = log(ZirconGrowth.zr_diffusivity(Th[i], XH2O))
        for j in 1:nelem
            kd_mat[i, j]  = log(kd[j])
            dif_mat[i, j] = log(dif[j])
        end
    end

    colors = _turbo_colors(nelem)

    # Helper: place a label on a line at a given x-position
    function _label_on_line!(ax, xv, yv, xpos, label; color=:black, fontsize=14)
        # Interpolate y at xpos
        idx = searchsortedfirst(xv, xpos)
        idx = clamp(idx, 2, length(xv))
        t = (xpos - xv[idx-1]) / (xv[idx] - xv[idx-1])
        ypos = yv[idx-1] + t * (yv[idx] - yv[idx-1])
        text!(ax, xpos, ypos; text=label, color=color, fontsize=fontsize,
              font=:bold, align=(:center, :bottom))
    end

    # Distribute label positions across the x range to reduce overlap
    xmin, xmax = extrema(xvals)
    n_labels = nelem + 1  # +1 for Zr
    label_xpos = range(xmin + 0.15 * (xmax - xmin),
                       xmax - 0.05 * (xmax - xmin), length=n_labels)

    ax1 = Axis(fig[1, 1]; xlabel="1000/T (K)", ylabel="ln Kd",
               title="Partition coefficients")
    lines!(ax1, xvals, zr_kd; color=:black, linewidth=3, linestyle=:dash)
    _label_on_line!(ax1, xvals, zr_kd, label_xpos[1], "Zr"; color=:black)
    for j in 1:nelem
        lines!(ax1, xvals, kd_mat[:, j]; color=colors[j], linewidth=2)
        _label_on_line!(ax1, xvals, kd_mat[:, j], label_xpos[j+1],
                        element_names[j]; color=colors[j])
    end

    ax2 = Axis(fig[1, 2]; xlabel="1000/T (K)", ylabel="ln D (cm²/s)",
               title="Diffusion coefficients")
    lines!(ax2, xvals, zr_dif; color=:black, linewidth=3, linestyle=:dash)
    _label_on_line!(ax2, xvals, zr_dif, label_xpos[1], "Zr"; color=:black)
    for j in 1:nelem
        lines!(ax2, xvals, dif_mat[:, j]; color=colors[j], linewidth=2)
        _label_on_line!(ax2, xvals, dif_mat[:, j], label_xpos[j+1],
                        element_names[j]; color=colors[j])
    end

    return fig
end

# ─── Figure 10 equivalent: element profiles for two timescales ────────────────

function ZirconGrowth.plot_element_profiles(short::SimulationResult{N},
                                            long::SimulationResult{N};
                                            elements::Vector{String} = ["U", "Y", "P", "Th", "Hf", "Sc"],
                                            short_label::String = "",
                                            long_label::String  = "",
                                            fig::Figure = Figure(size=(1200, 600))) where {N}
    nr_s = _normalize_radius(short.zircon_radius_um)
    nr_l = _normalize_radius(long.zircon_radius_um)

    zc_s = short.interface_melt .* short.partition_history
    zc_l = long.interface_melt  .* long.partition_history

    for (idx, elem) in enumerate(elements)
        row = (idx - 1) ÷ 3 + 1
        col = (idx - 1) % 3 + 1
        ie_s = _find_element(short.element_names, elem)
        ie_l = _find_element(long.element_names, elem)
        ie_s === nothing && continue

        ax = Axis(fig[row, col]; xlabel="Normalised radius",
                  ylabel="$elem (ppm)", title=elem)
        lines!(ax, nr_s, zc_s[:, ie_s]; linewidth=2, label=short_label)
        lines!(ax, nr_l, zc_l[:, ie_l]; linewidth=2, label=long_label)
        axislegend(ax; position=:rb, framevisible=false)
    end
    return fig
end

# ─── Figure 11 equivalent: element ratios ────────────────────────────────────

function ZirconGrowth.plot_element_ratios(short::SimulationResult{N},
                                          long::SimulationResult{N};
                                          ratios::Vector{Tuple{String,String}} = [
                                              ("Y","Th"), ("Y","Hf"), ("Hf","P"),
                                              ("Th","U"), ("Y","Sc"), ("Sc","Nb")],
                                          short_label::String = "",
                                          long_label::String  = "",
                                          fig::Figure = Figure(size=(1200, 600))) where {N}
    nr_s = _normalize_radius(short.zircon_radius_um)
    nr_l = _normalize_radius(long.zircon_radius_um)

    zc_s = short.interface_melt .* short.partition_history
    zc_l = long.interface_melt  .* long.partition_history

    for (idx, (num, den)) in enumerate(ratios)
        row = (idx - 1) ÷ 3 + 1
        col = (idx - 1) % 3 + 1

        in_s = _find_element(short.element_names, num)
        id_s = _find_element(short.element_names, den)
        in_l = _find_element(long.element_names,  num)
        id_l = _find_element(long.element_names,  den)
        (in_s === nothing || id_s === nothing) && continue

        r_s = _safe_ratio(zc_s[:, in_s], zc_s[:, id_s])
        r_l = _safe_ratio(zc_l[:, in_l], zc_l[:, id_l])

        ax = Axis(fig[row, col]; xlabel="Normalised radius",
                  ylabel="$num/$den", title="$num/$den")
        lines!(ax, nr_s, r_s; linewidth=2, label=short_label)
        lines!(ax, nr_l, r_l; linewidth=2, label=long_label)
        axislegend(ax; position=:rb, framevisible=false)
    end
    return fig
end

# ─── Growth-diagnostic plots ─────────────────────────────────────────────────

function ZirconGrowth.plot_growth_diagnostics(res::SimulationResult;
                                              fig::Figure = Figure(size=(900, 700)))
    ax1 = Axis(fig[1, 1]; xlabel="r (µm)", ylabel="C (ppm)", title="Zr melt profile")
    lines!(ax1, res.radius_um, res.concentrations[:, 1]; linewidth=2)

    ax2 = Axis(fig[1, 2]; xlabel="Time (yr)", ylabel="Zr radius (µm)",
               title="Crystal radius")
    lines!(ax2, res.time_years, res.zircon_radius_um; linewidth=2)

    ax3 = Axis(fig[2, 1]; xlabel="Time (yr)", ylabel="µm yr⁻¹",
               title="Growth & contraction rates")
    lines!(ax3, res.time_years, res.growth_rate; linewidth=2, label="V (growth)")
    lines!(ax3, res.time_years, .-res.matrix_rate; linewidth=2, label="-W (contraction)")
    axislegend(ax3; position=:rt, framevisible=false)

    ax4 = Axis(fig[2, 2]; xlabel="Zr radius (µm)", ylabel="C/Cmax",
               title="Normalised rim concentrations")
    ntr = size(res.rim_conc, 2)
    colors = _turbo_colors(ntr)
    for j in 1:ntr
        cmax = maximum(view(res.rim_conc, 2:size(res.rim_conc, 1), j))
        cmax == 0.0 && continue
        lines!(ax4, res.zircon_radius_um[2:end],
               res.rim_conc[2:end, j] ./ cmax;
               linewidth=1.2, color=colors[j], label=res.element_names[j])
    end
    axislegend(ax4; position=:lt, framevisible=false)

    return fig
end

# ─── Helpers ──────────────────────────────────────────────────────────────────

function _normalize_radius(r::AbstractVector{Float64})
    rng = r[end] - r[1]
    if abs(rng) < eps(maximum(abs, r) + 1.0)
        return zeros(length(r))
    end
    return (r .- r[1]) ./ rng
end

function _find_element(names::NTuple{N,String}, target::String) where {N}
    tl = lowercase(target)
    for i in 1:N
        lowercase(names[i]) == tl && return i
    end
    return nothing
end

function _safe_ratio(num::AbstractVector, den::AbstractVector)
    out = similar(num)
    for i in eachindex(num, den)
        d = den[i]
        out[i] = abs(d) < eps(maximum(abs, den) + 1.0) ? NaN : num[i] / d
    end
    return out
end

end # module
