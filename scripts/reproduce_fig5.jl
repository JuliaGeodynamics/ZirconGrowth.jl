#!/usr/bin/env julia
"""
    reproduce_fig5.jl

Reproduce Figure 5 from:
  Bindeman & Melnik (2026) Earth Planet. Sci. Lett. 673, 119723.

Figure 5 is a 2×3 panel showing trace-element concentration profiles
across zircon (normalised radius: core = 0, rim = 1) for six selected
elements that span a wide range of Kd and D values:

  Row 1: U  (highly compatible)    | Y  (very highly compatible) | P  (moderately compatible)
  Row 2: Th (highly compatible)    | Hf (very highly compatible) | Sc (highly incompatible)

Each panel compares fast (500 yr, ~43 µm) vs slow (50 000 yr, ~145 µm)
crystallisation using default parameters: 940→680 °C, 2 wt% H₂O, M = 1.3.

Run from the package root:
    julia --project scripts/reproduce_fig5.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ZirconGrowth
using GLMakie

function run_fig5()
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

    # ── Elements to plot (same order as in the paper) ────────────────────
    elements = ["U", "Y", "P", "Th", "Hf", "Sc"]
    titles   = [
        "U (D: small, Kd high)",
        "Y (D: medium, Kd high)",
        "P (D: small, Kd small)",
        "Th (D: small, Kd medium)",
        "Hf (D: high, Kd highest)",
        "Sc (D: small, Kd very small)",
    ]
    ylimits = Dict(
        "U"  => (0, 600),
        "Y"  => (0, 12000),
        "P"  => (0, 14000),
        "Th" => (100, 700),
        "Hf" => (0, 8000),
        "Sc" => (0, 40),
    )

    # Map element names to column indices
    function elem_idx(names, el)
        tl = lowercase(el)
        for i in eachindex(names)
            lowercase(names[i]) == tl && return i
        end
        error("Element $el not found")
    end

    # ── Figure ───────────────────────────────────────────────────────────
    fig = Figure(size=(1200, 600))

    for (idx, (el, ttl)) in enumerate(zip(elements, titles))
        row = (idx - 1) ÷ 3 + 1
        col = (idx - 1) % 3 + 1

        ie = elem_idx(res_fast.element_names, el)

        ax = Axis(fig[row, col];
                  xlabel = row == 2 ? "Normalised radius" : "",
                  ylabel = "$el (ppm)",
                  title  = ttl,
                  limits = (nothing, ylimits[el]))

        lines!(ax, nr_fast, zc_fast[:, ie]; linewidth=2, label="500 yr")
        lines!(ax, nr_slow, zc_slow[:, ie]; linewidth=2, label="50 000 yr")
        axislegend(ax; position=:lt, framevisible=false)
    end

    # Print final radii for comparison with paper values (43 µm, 145 µm)
    @info "Final zircon radii (paper: ≈43 µm @ 500 yr, ≈145 µm @ 50 kyr)" r_fast=res_fast.zircon_radius_um[end] r_slow=res_slow.zircon_radius_um[end]

    return fig
end

fig = run_fig5()
display(fig)
@info "Figure 5 displayed."
