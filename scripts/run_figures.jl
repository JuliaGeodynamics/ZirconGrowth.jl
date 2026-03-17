#!/usr/bin/env julia
"""
    run_figures.jl

Reproduce the figures from the MATLAB code using ZirconGrowth.jl and GLMakie.
Run from the package root: `julia --project scripts/run_figures.jl`
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ZirconGrowth
using GLMakie  # triggers ZirconGrowthGLMakieExt

# ─── Run simulations ──────────────────────────────────────────────────────────
@info "Running 500-year simulation…"
res_short = simulate_zircon_growth(500)

@info "Running 50 000-year simulation…"
res_long  = simulate_zircon_growth(50_000)

# ─── Figure 1: Kd and D vs 1000/T ────────────────────────────────────────────
Td = 940.0 + 273.0
T0 = 680.0 + 273.0
nt_plot = 500
Th_plot = [ZirconGrowth.temperature_history(Td, T0, 0, 0.0, (i-1) * 1.0 / (nt_plot-1), 1.0)
           for i in 1:nt_plot]

fig1 = plot_trace_properties(Th_plot, 2.0, default_element_data().names)
display(fig1)

# ─── Figure 2: Element profiles (500 yr vs 50 000 yr) ────────────────────────
fig2 = plot_element_profiles(res_short, res_long;
                             elements=["U", "Y", "P", "Th", "Hf", "Sc"],
                             short_label="500 yr", long_label="50 000 yr")
display(fig2)

# ─── Figure 3: Element ratios ────────────────────────────────────────────────
fig3 = plot_element_ratios(res_short, res_long;
                           short_label="500 yr", long_label="50 000 yr")
display(fig3)

# ─── Figure 4: Growth diagnostics (500 yr) ──────────────────────────────────
fig4 = plot_growth_diagnostics(res_short)
display(fig4)

@info "All figures displayed."
