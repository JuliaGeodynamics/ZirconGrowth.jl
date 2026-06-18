#!/usr/bin/env julia
"""
    custom_cooling_path.jl

Example: run ZirconGrowth.jl with a user-supplied, irregularly-spaced
cooling history given in °C and Myr.

The workflow is:
  1. Define (time_Myr, T_Celsius) data points (may be irregularly spaced).
  2. Convert to (time_yr, T_K).
  3. Interpolate onto the uniform internal time grid used by the solver.
  4. Pass the result as `Th_override` to `simulate_zircon_growth!`.

No external packages are required: interpolation is done with a simple
piecewise-linear helper below.

Run from the package root:
    julia --project scripts/custom_cooling_path.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ZirconGrowth
using Printf

# ─── User data: irregularly-spaced (time [Myr], T [°C]) ───────────────────────
#
# Two-stage cooling with a mid-path reheating episode (e.g. magma recharge):
#   0 – 0.03 Myr  : rapid post-emplacement cooling  (940 → 810 °C)
#   0.03 – 0.05 Myr: reheating (magma recharge)      (810 → 860 °C)
#   0.05 – 0.20 Myr: final cooling to solidus         (860 → 680 °C)
#
# Data points are deliberately irregularly spaced to illustrate the workflow.
# The temperature path may be non-monotonic; the solver handles it correctly.

time_Myr = Float64[0.000, 0.005, 0.015, 0.030, 0.040, 0.050, 0.075, 0.150, 0.200]
T_C      = Float64[940,   920,   870,   810,   840,   860,   790,   720,   680]

# ─── Simulation parameters ────────────────────────────────────────────────────

p = GrowthParams(
    tfinal_years = time_Myr[end] * 1e6,
    nt           = 5000,
    nx           = 500,
    Td           = T_C[1]   + 273.15,
    T0           = T_C[end] + 273.15,
    XH2O         = 2.0,
    Mfactor      = 1.3,
)

# ─── Run simulation ───────────────────────────────────────────────────────────

result = simulate_from_cooling_path(time_Myr, T_C; params=p)

# ─── Inspect results ──────────────────────────────────────────────────────────

@info "Final zircon radius" radius_um=result.zircon_radius_um[end]
@info "Time span"           from_yr=result.time_years[1] to_yr=result.time_years[end]

# Print rim concentrations for selected elements
println("\nRim concentrations at final time step (ppm):")
for (i, el) in enumerate(result.element_names)
    @printf "  %-4s  %8.2f ppm\n" el result.rim_conc[end, i]
end
