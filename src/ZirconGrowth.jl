"""
    ZirconGrowth

Numerical model of zircon crystal growth in a cooling magma, tracking the
evolution of Zr saturation and trace-element partitioning.  The model solves
a moving-boundary diffusion problem on a 1-D spherical grid (radius) using
a Crank–Nicolson-like implicit scheme with a tridiagonal (Thomas) solver.

# Physical picture
A zircon crystal of radius `s` sits at the centre of a spherical melt cell
of radius `R`.  As the melt cools, Zr saturation drops and the crystal
grows; at the same time, 11 trace elements (Hf, Y, U, Th, Sm, Dy, Yb, P,
Nb, Sc, Li) diffuse through the melt and are incorporated into the crystal
rim according to temperature-dependent partition coefficients ``K_d(T)``.

# Key references
* Boehnke et al. (2013) Chem. Geol. 351, 324 — Zr saturation model
* Watson (1996) — Zr diffusion parameterisation
"""
module ZirconGrowth

using DelimitedFiles

export ElementData, GrowthParams, SimulationResult, SimulationWorkspace
export load_element_data, default_element_data
export zr_saturation, zr_diffusivity, trace_diffusivities, partition_coefficients
export temperature_history, plagioclase_fraction, plagioclase_fraction_poly
export simulate_zircon_growth, simulate_zircon_growth!, compute_profiles
export simulate_from_cooling_path
export lerp, lerp_vec
export plot_trace_properties, plot_element_profiles, plot_element_ratios,
       plot_growth_diagnostics

# ─── plotting stubs (methods provided by ext/ZirconGrowthGLMakieExt.jl) ───────
function plot_trace_properties end
function plot_element_profiles end
function plot_element_ratios end
function plot_growth_diagnostics end

# ─── data directory ────────────────────────────────────────────────────────────
const DATA_DIR = joinpath(@__DIR__, "..", "data")

# ─── source files ──────────────────────────────────────────────────────────────
include("types.jl")
include("physics.jl")
include("solver.jl")
include("simulation.jl")

end # module
