# ZirconGrowth.jl

[![CI](https://github.com/JuliaGeodynamics/ZirconGrowth.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaGeodynamics/ZirconGrowth.jl/actions/workflows/CI.yml)

This Julia package simulates zircon growth from a melt along with trace element evolution. It is a Julia translation of the MATLAB codes accompanying:

> Bindeman, I., Melnik, O., 2026. Non-equilibrium behavior of trace elements during zircon crystallization from the melt. *Earth and Planetary Science Letter*s* 673, 119723. https://doi.org/10.1016/j.epsl.2025.119723

Original MATLAB implementation by **Oleg Melnik**. This package is a line-for-line translation into Julia, validated against the MATLAB output to 8+ significant figures.

Numerical model of zircon crystal growth in a cooling silicate melt, tracking Zr saturation and 11 trace-element (Hf, Y, U, Th, Sm, Dy, Yb, P, Nb, Sc, Li) concentration profiles.

## Quick start

```julia
using ZirconGrowth

# Run a 500-year cooling simulation (default parameters)
result = simulate_zircon_growth(500)

# Inspect output
result.zircon_radius_um   # zircon radius history (µm)
result.rim_conc           # nt × 11 matrix of rim concentrations (ppm)
result.time_years         # time array (yr)
result.element_names      # ("Hf", "Y", "U", …)
```

### Custom parameters

```julia
p = GrowthParams(
    tfinal_years = 1000.0,
    nt   = 1000,          # time steps
    nx   = 300,           # spatial grid points
    Td   = 950.0 + 273,   # initial T (K)
    T0   = 700.0 + 273,   # final T (K)
    XH2O = 3.0,           # water content (wt%)
)
result = simulate_zircon_growth(1000; params=p)
```

### Allocation-free in-place interface

Pre-allocate a workspace and reuse it for zero-allocation repeated runs:

```julia
using ZirconGrowth

p  = GrowthParams(tfinal_years=500.0, nt=500, nx=200)
ed = default_element_data()
ws = SimulationWorkspace(p, ed)

# First run — zero heap allocations in the hot loop
result = simulate_zircon_growth!(ws, 500.0, p, ed)

# Subsequent runs reuse the same workspace
result = simulate_zircon_growth!(ws, 500.0, p, ed)  # still 0 allocations
```

### Plotting (GLMakie extension)

Plotting functions are provided via a package extension that loads
automatically when GLMakie is imported — no `include` needed:

```julia
using ZirconGrowth
using GLMakie   # triggers the plotting extension

res_short = simulate_zircon_growth(500)
res_long  = simulate_zircon_growth(50_000)

fig = plot_element_profiles(res_short, res_long;
        short_label="500 yr", long_label="50 000 yr")
```

## Package structure

```
ZirconGrowth.jl/
├── Project.toml
├── src/
│   ├── ZirconGrowth.jl   # module definition
│   ├── types.jl           # ElementData, GrowthParams, SimulationResult
│   ├── physics.jl         # T-dependent Zr saturation, diffusion, Kd
│   ├── solver.jl          # Thomas tridiagonal solver, coefficient assembly
│   └── simulation.jl      # main time-stepping loop
├── test/
│   └── runtests.jl        # 56 tests incl. MATLAB reference validation
├── ext/
│   └── ZirconGrowthGLMakieExt.jl  # GLMakie plotting extension
├── data/
│   └── element_data.csv   # default trace-element concentrations
└── scripts/
    ├── plot_recipes.jl    # reusable plotting helpers
    ├── reproduce_fig3.jl  # reproduce Figure 3 from the paper
    └── run_figures.jl     # reproduce all paper figures
```

## Design goals

- **Allocation-free inner loop** — all workspace vectors are pre-allocated;
  the main time loop performs zero heap allocations.
- **Validated against MATLAB** — the test suite compares final zircon radius
  and rim concentrations against the reference MATLAB output to 8+ significant
  figures.
- **Clean, documented API** — every public function and type has a docstring.

## Running tests

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## References

- Bindeman, I., Melnik, O., 2026. Non-equilibrium behavior of trace elements during zircon crystallization from the melt. *Earth and Planetary Science Letters* 673, 119723. https://doi.org/10.1016/j.epsl.2025.119723
- Boehnke et al. (2013) *Chem. Geol.* **351**, 324 — Zr saturation model
- Watson (1996) — Zr diffusion parameterisation
