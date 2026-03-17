#!/usr/bin/env julia
"""
    reproduce_fig2.jl

Reproduce Figure 2 from:
  Bindeman & Melnik (2026) Earth Planet. Sci. Lett. 673, 119723.

Figure 2 shows temperature-dependent partition coefficients (A) and
diffusion coefficients (B) as ln(Kd) and ln(D) vs 1000/T for the 11
trace elements (Hf, Y, U, Th, Sm, Dy, Yb, P, Nb, Sc, Li) plus Zr
(dashed). Temperature range: 940–680 °C (1213–953 K), 2 wt% H₂O.

Run from the package root:
    julia --project scripts/reproduce_fig2.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ZirconGrowth
using GLMakie

function run_fig2()
    Td = 940.0 + 273.0
    T0 = 680.0 + 273.0
    nt = 500

    # Temperature array spanning the cooling range
    Th = [ZirconGrowth.temperature_history(Td, T0, 0, 0.0,
            (i - 1) / (nt - 1), 1.0) for i in 1:nt]

    fig = plot_trace_properties(Th, 2.0, default_element_data().names;
                                fig=Figure(size=(1200, 500)))

    # Match paper axis limits
    ax_a = content(fig[1, 1])
    ax_b = content(fig[1, 2])
    ylims!(ax_a, -4, 10)
    ylims!(ax_b, -37, -22)

    return fig
end

fig = run_fig2()
display(fig)
@info "Figure 2 displayed."
