#!/usr/bin/env julia
"""
    reproduce_fig6.jl

Reproduce Figure 6 from:
  Bindeman & Melnik (2026) Earth Planet. Sci. Lett. 673, 119723.

Figure 6 shows trace-element ratios across zircon (normalised radius:
core = 0, rim = 1), comparing fast (500 yr) vs slow (50 000 yr)
crystallisation in a 3×3 grid:

  Row 1: Y/Th  | Y/Hf  | Hf/P
  Row 2: Th/U  | Y/Sc  | Sc/Nb
  Row 3: Zr/Hf | (δ94/90Zr — not implemented, see note)

The paper also shows δ94/90Zr (from Bindeman & Melnik 2022), which
requires an isotopic fractionation model not included in this package.

Default parameters: 940→680 °C, 2 wt% H₂O, M = 1.3.

Run from the package root:
    julia --project scripts/reproduce_fig6.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ZirconGrowth
using GLMakie

function run_fig6()
    # ── Simulations ──────────────────────────────────────────────────────
    p_fast = GrowthParams(tfinal_years=500.0,    nt=500,  nx=500)
    p_slow = GrowthParams(tfinal_years=50_000.0, nt=5000, nx=500)

    res_fast = simulate_zircon_growth(500;     params=p_fast)
    res_slow = simulate_zircon_growth(50_000;  params=p_slow)

    # ── Derived: zircon concentration = interface_melt × Kd ──────────────
    zc_fast = res_fast.interface_melt .* res_fast.partition_history
    zc_slow = res_slow.interface_melt .* res_slow.partition_history

    # Normalise radius (core = 0, rim = 1)
    function norm_radius(r)
        rng = r[end] - r[1]
        return abs(rng) < eps(maximum(abs, r) + 1.0) ? zeros(length(r)) :
               (r .- r[1]) ./ rng
    end
    nr_fast = norm_radius(res_fast.zircon_radius_um)
    nr_slow = norm_radius(res_slow.zircon_radius_um)

    # Map element names to column indices
    function elem_idx(names, el)
        tl = lowercase(el)
        for i in eachindex(names)
            lowercase(names[i]) == tl && return i
        end
        error("Element $el not found")
    end

    # Safe ratio avoiding division by zero
    function safe_ratio(num, den)
        out = similar(num)
        thresh = eps(maximum(abs, den) + 1.0)
        for i in eachindex(num, den)
            out[i] = abs(den[i]) < thresh ? NaN : num[i] / den[i]
        end
        return out
    end

    # ── Ratios to plot (same as MATLAB main.m / paper Fig. 6) ────────────
    ratios = [
        ("Y",  "Th"),
        ("Y",  "Hf"),
        ("Hf", "P"),
        ("Th", "U"),
        ("Y",  "Sc"),
        ("Sc", "Nb"),
    ]

    # ── Figure ───────────────────────────────────────────────────────────
    fig = Figure(size=(1200, 900))

    for (idx, (num, den)) in enumerate(ratios)
        row = (idx - 1) ÷ 3 + 1
        col = (idx - 1) % 3 + 1

        in_f = elem_idx(res_fast.element_names, num)
        id_f = elem_idx(res_fast.element_names, den)
        in_s = elem_idx(res_slow.element_names, num)
        id_s = elem_idx(res_slow.element_names, den)

        r_fast = safe_ratio(zc_fast[:, in_f], zc_fast[:, id_f])
        r_slow = safe_ratio(zc_slow[:, in_s], zc_slow[:, id_s])

        ax = Axis(fig[row, col];
                  xlabel = row == 3 ? "Normalised radius" : "",
                  ylabel = "$num/$den",
                  title  = "$num/$den")

        lines!(ax, nr_fast, r_fast; linewidth=2, label="500 yr")
        lines!(ax, nr_slow, r_slow; linewidth=2, label="50 000 yr")
        axislegend(ax; position=:lt, framevisible=false)
    end

    # ── Panel 7: Zr/Hf ──────────────────────────────────────────────────
    # Zr in zircon is the constant Czirc; Hf in zircon = interface_melt_Hf × Kd_Hf
    Czirc = 490_000.0
    ihf_f = elem_idx(res_fast.element_names, "Hf")
    ihf_s = elem_idx(res_slow.element_names, "Hf")
    zrhf_fast = safe_ratio(fill(Czirc, size(zc_fast, 1)), zc_fast[:, ihf_f])
    zrhf_slow = safe_ratio(fill(Czirc, size(zc_slow, 1)), zc_slow[:, ihf_s])

    ax_zrhf = Axis(fig[3, 1];
                   xlabel = "Normalised radius",
                   ylabel = "Zr/Hf",
                   title  = "Zr/Hf")
    lines!(ax_zrhf, nr_fast, zrhf_fast; linewidth=2, label="500 yr")
    lines!(ax_zrhf, nr_slow, zrhf_slow; linewidth=2, label="50 000 yr")
    axislegend(ax_zrhf; position=:lt, framevisible=false)

    # ── Panel 8: δ94/90Zr (placeholder) ─────────────────────────────────
    # The isotopic fractionation model from Bindeman & Melnik (2022) is not
    # included in this package. Show an annotation instead.
    ax_d = Axis(fig[3, 2];
                xlabel = "Normalised radius",
                ylabel = "δ⁹⁴/⁹⁰Zr (‰)",
                title  = "δ⁹⁴/⁹⁰Zr")
    text!(ax_d, 0.5, 0.5; text="Not implemented\n(see Bindeman &\nMelnik, 2022)",
          align=(:center, :center), fontsize=14, space=:relative)

    return fig
end

fig = run_fig6()
display(fig)
@info "Figure 6 displayed."
