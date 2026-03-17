#!/usr/bin/env julia
"""
    reproduce_fig3.jl

Reproduce Figure 3 from:
  Bindeman & Melnik (2026) Earth Planet. Sci. Lett. 673, 119723.

Figure 3 is a 2×2 panel comparing 500-year (fast) and 50,000-year (slow)
zircon crystallisation, each with and without major-mineral fractionation:
  (A) Zircon radius vs time              [upper-left]
  (B) Zircon volume vs time               [lower-left]
  (C) Radial growth rate vs zircon radius [upper-right]
  (D) Radial growth rate vs time          [lower-right]

Solid lines:  without major minerals
Dashed lines: with major minerals (Kd(Zr) = 1)

Paper parameters: cooling from 930 °C to 680 °C, Rc = 0.2 cm (2 mm),
s₀ = 1 µm, 2 wt% H₂O, M = 1.3.

Run from the package root:
    julia --project scripts/reproduce_fig3.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ZirconGrowth
using GLMakie

function run_fig3()
    # Paper parameters: 930 °C to 680 °C (Td = 1203 K, T0 = 953 K)
    Td = 930.0 + 273.0
    T0 = 680.0 + 273.0

    # -- Run four simulations -------------------------------------------------
    # Without major minerals (solid lines)
    p_fast    = GrowthParams(tfinal_years=500.0,    nt=500,  nx=500, Td=Td, T0=T0)
    p_slow    = GrowthParams(tfinal_years=50_000.0, nt=5000, nx=500, Td=Td, T0=T0)
    # With major minerals, Kd(Zr) = 1 (dashed lines)
    p_fast_mm = GrowthParams(tfinal_years=500.0,    nt=500,  nx=500, Td=Td, T0=T0,
                             major_minerals=true, kpl=1.0)
    p_slow_mm = GrowthParams(tfinal_years=50_000.0, nt=5000, nx=500, Td=Td, T0=T0,
                             major_minerals=true, kpl=1.0)

    res_fast    = simulate_zircon_growth(500;     params=p_fast)
    res_slow    = simulate_zircon_growth(50_000;  params=p_slow)
    res_fast_mm = simulate_zircon_growth(500;     params=p_fast_mm)
    res_slow_mm = simulate_zircon_growth(50_000;  params=p_slow_mm)

    # -- Helper: derived quantities for one result ----------------------------
    function derived(res)
        nt = length(res.time_years)
        # Radial growth rate: µm/yr → cm/s
        um_yr_to_cm_s = 1e-4 / (365.25 * 24 * 3600)
        dr_dt = max.(res.growth_rate .* um_yr_to_cm_s, 1e-30)  # clamp for log

        # Volume: V = (4/3)πr³  [µm³]
        vol = (4π / 3) .* res.zircon_radius_um.^3

        return (; nt, dr_dt, vol)
    end

    df  = derived(res_fast)
    ds  = derived(res_slow)
    dfm = derived(res_fast_mm)
    dsm = derived(res_slow_mm)

    # Common log-scale time ticks
    time_ticks = [1, 10, 100, 1_000, 10_000]
    time_ticklabels = ["1", "10", "100", "1,000", "10,000"]

    # -- Figure ---------------------------------------------------------------
    fig = Figure(size=(1000, 800))

    # (A) Zircon radius vs time [upper-left]
    ax_a = Axis(fig[1, 1]; xlabel="Time (years)",
                ylabel="Zircon radius (µm)", title="A",
                xscale=log10,
                xticks=(time_ticks, time_ticklabels))
    lines!(ax_a, res_fast.time_years[2:df.nt], res_fast.zircon_radius_um[2:df.nt];
           linewidth=2, label="500 yr")
    lines!(ax_a, res_slow.time_years[2:ds.nt], res_slow.zircon_radius_um[2:ds.nt];
           linewidth=2, label="50,000 yr")
    lines!(ax_a, res_fast_mm.time_years[2:dfm.nt], res_fast_mm.zircon_radius_um[2:dfm.nt];
           linewidth=2, linestyle=:dash, label="500 yr + minerals")
    lines!(ax_a, res_slow_mm.time_years[2:dsm.nt], res_slow_mm.zircon_radius_um[2:dsm.nt];
           linewidth=2, linestyle=:dash, label="50,000 yr + minerals")
    axislegend(ax_a; position=:lt, framevisible=false)

    # (C) Radial growth rate vs zircon radius [upper-right]
    ax_c = Axis(fig[1, 2]; xlabel="Zircon radius (µm)",
                ylabel="Growth rate (cm/s)", title="C", yscale=log10,
                limits=(nothing, (1e-16, 1e-11)))
    lines!(ax_c, res_fast.zircon_radius_um[2:df.nt], df.dr_dt[2:df.nt];
           linewidth=2, label="500 yr")
    lines!(ax_c, res_slow.zircon_radius_um[2:ds.nt], ds.dr_dt[2:ds.nt];
           linewidth=2, label="50,000 yr")
    lines!(ax_c, res_fast_mm.zircon_radius_um[2:dfm.nt], dfm.dr_dt[2:dfm.nt];
           linewidth=2, linestyle=:dash, label="500 yr + minerals")
    lines!(ax_c, res_slow_mm.zircon_radius_um[2:dsm.nt], dsm.dr_dt[2:dsm.nt];
           linewidth=2, linestyle=:dash, label="50,000 yr + minerals")
    axislegend(ax_c; position=:rt, framevisible=false)

    # (B) Zircon volume vs time [lower-left]
    ax_b = Axis(fig[2, 1]; xlabel="Time (years)",
                ylabel="Volume (µm³)", title="B",
                xscale=log10, yscale=log10,
                xticks=(time_ticks, time_ticklabels))
    lines!(ax_b, res_fast.time_years[2:df.nt], df.vol[2:df.nt];
           linewidth=2, label="500 yr")
    lines!(ax_b, res_slow.time_years[2:ds.nt], ds.vol[2:ds.nt];
           linewidth=2, label="50,000 yr")
    lines!(ax_b, res_fast_mm.time_years[2:dfm.nt], dfm.vol[2:dfm.nt];
           linewidth=2, linestyle=:dash, label="500 yr + minerals")
    lines!(ax_b, res_slow_mm.time_years[2:dsm.nt], dsm.vol[2:dsm.nt];
           linewidth=2, linestyle=:dash, label="50,000 yr + minerals")
    axislegend(ax_b; position=:rb, framevisible=false)

    # (D) Radial growth rate vs time [lower-right]
    ax_d = Axis(fig[2, 2]; xlabel="Time (years)",
                ylabel="Growth rate (cm/s)", title="D",
                xscale=log10, yscale=log10,
                xticks=(time_ticks, time_ticklabels),
                limits=(nothing, (1e-16, 1e-11)))
    lines!(ax_d, res_fast.time_years[2:df.nt], df.dr_dt[2:df.nt];
           linewidth=2, label="500 yr")
    lines!(ax_d, res_slow.time_years[2:ds.nt], ds.dr_dt[2:ds.nt];
           linewidth=2, label="50,000 yr")
    lines!(ax_d, res_fast_mm.time_years[2:dfm.nt], dfm.dr_dt[2:dfm.nt];
           linewidth=2, linestyle=:dash, label="500 yr + minerals")
    lines!(ax_d, res_slow_mm.time_years[2:dsm.nt], dsm.dr_dt[2:dsm.nt];
           linewidth=2, linestyle=:dash, label="50,000 yr + minerals")
    axislegend(ax_d; position=:lb, framevisible=false)

    # Print final radii for comparison with paper values
    @info "Final radii (paper: 43 µm @ 500 yr, 145 µm @ 50 kyr)" r_fast=res_fast.zircon_radius_um[end] r_slow=res_slow.zircon_radius_um[end] r_fast_mm=res_fast_mm.zircon_radius_um[end] r_slow_mm=res_slow_mm.zircon_radius_um[end]

    return fig
end

fig = run_fig3()
display(fig)
@info "Figure 3 displayed."
