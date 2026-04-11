# ─── Element data ──────────────────────────────────────────────────────────────

"""
    ElementData

Immutable container for the names and initial melt concentrations (ppm)
of the 11 trace elements tracked by the model.

# Fields
- `names  :: NTuple{N,String}` — element symbols (e.g. `"Hf"`, `"Y"`, …)
- `concentrations :: NTuple{N,Float64}` — initial melt concentrations in ppm

The element order is fixed and must match the ordering assumed by
[`partition_coefficients`](@ref) and [`trace_diffusivities`](@ref):
Hf, Y, U, Th, Sm, Dy, Yb, P, Nb, Sc, Li.
"""
struct ElementData{N}
    names::NTuple{N,String}
    concentrations::NTuple{N,Float64}
end

Base.length(::ElementData{N}) where {N} = N

"""
    load_element_data(path::String) -> ElementData

Read an `element_data.csv` file (two columns: `Element`, `Concentration`)
and return an [`ElementData`](@ref) instance.
"""
function load_element_data(path::String)
    lines = readlines(path)
    names  = String[]
    concs  = Float64[]
    for line in lines[2:end]           # skip header
        isempty(strip(line)) && continue
        parts = split(line, ',')
        push!(names, strip(String(parts[1])))
        push!(concs, parse(Float64, parts[2]))
    end
    N = length(names)
    return ElementData{N}(NTuple{N,String}(names), NTuple{N,Float64}(concs))
end

"""
    default_element_data() -> ElementData{11}

Load the built-in element data shipped with the package.
"""
function default_element_data()
    load_element_data(joinpath(DATA_DIR, "element_data.csv"))
end

# ─── Growth parameters ────────────────────────────────────────────────────────

"""
    GrowthParams

All physical and numerical parameters needed to run a simulation.

# Numerical
- `nt::Int` — number of time steps
- `nx::Int` — number of spatial grid points

# Physical
- `Rc::Float64`    — melt-cell radius (cm), default 0.2
- `s0::Float64`    — initial zircon radius (cm), default 1e-4
- `Td::Float64`    — initial (liquidus) temperature (K), default 1213
- `T0::Float64`    — final (solidus) temperature (K), default 953
- `Czirc::Float64` — Zr concentration in zircon (ppm), default 490 000
- `kpl::Float64`   — bulk partition coefficient for Zr in major phases, default 0.2
- `XH2O::Float64`  — water content (wt%), default 2.0
- `Mfactor::Float64` — M-factor for Boehnke saturation model, default 1.3

# Temperature history
- `n_oscillations::Int` — number of T-oscillation cycles, default 0
- `amplitude::Float64`  — T-oscillation amplitude (K), default 0.0
"""
struct GrowthParams
    # numerical
    nt::Int
    nx::Int
    # physical
    Rc::Float64
    s0::Float64
    Td::Float64
    T0::Float64
    Czirc::Float64
    kpl::Float64
    XH2O::Float64
    Mfactor::Float64
    # temperature history
    n_oscillations::Int
    amplitude::Float64
    # major-mineral fractionation
    major_minerals::Bool
end

"""
    GrowthParams(; tfinal_years, kwargs...)

Convenience constructor.  `nt` defaults to `clamp(tfinal_years, 500, 5000)` when
not provided, and `nx = 500`.
"""
function GrowthParams(;
    tfinal_years::Float64,
    nt::Int = clamp(round(Int, tfinal_years), 500, 5000),
    nx::Int = 500,
    Rc::Float64 = 0.2,
    s0::Float64 = 1e-4,
    Td::Float64 = 940.0 + 273.0,
    T0::Float64 = 680.0 + 273.0,
    Czirc::Float64 = 490_000.0,
    kpl::Float64 = 0.2,
    XH2O::Float64 = 2.0,
    Mfactor::Float64 = 1.3,
    n_oscillations::Int = 0,
    amplitude::Float64 = 0.0,
    major_minerals::Bool = false,
)
    return GrowthParams(nt, nx, Rc, s0, Td, T0, Czirc, kpl, XH2O, Mfactor,
                        n_oscillations, amplitude, major_minerals)
end

"""
    GrowthParams(time_Myr, T_C; kwargs...)

Construct a `GrowthParams` from a cooling-path description.  `tfinal_years`,
`Td`, and `T0` are derived automatically from the input vectors:

- `tfinal_years = time_Myr[end] * 1e6`
- `Td  = T_C[1]   + 273.15`   (starting temperature)
- `T0  = T_C[end] + 273.15`   (ending temperature)

Any keyword accepted by `GrowthParams(; tfinal_years, ...)` can be passed to
override the derived values or set other parameters.

# Example
```julia
time_Myr = [0.0, 0.05, 0.1, 0.2]
T_C      = [940.0, 860.0, 780.0, 680.0]
p = GrowthParams(time_Myr, T_C; nx=300, nt=2000)
```
"""
function GrowthParams(time_Myr::AbstractVector{<:Real}, T_C::AbstractVector{<:Real}; kwargs...)
    return GrowthParams(;
        tfinal_years = Float64(time_Myr[end]) * 1e6,
        Td           = Float64(T_C[1])   + 273.15,
        T0           = Float64(T_C[end]) + 273.15,
        kwargs...
    )
end

# ─── Simulation result ────────────────────────────────────────────────────────

"""
    SimulationResult

Output container returned by [`simulate`](@ref).

# Time-series fields (length `nt`)
- `time_years`       — physical time (years)
- `zircon_radius_um` — zircon radius (µm)
- `growth_rate`      — crystal growth velocity (µm yr⁻¹)
- `matrix_rate`      — melt-cell contraction rate (µm yr⁻¹)
- `rim_conc`         — `nt × ntr` matrix of rim concentrations (ppm)
- `partition_history` — `nt × ntr` matrix of ``K_d`` at each step
- `interface_melt`   — `nt × ntr` matrix of melt concentration at the interface

# Spatial profile fields (length `nx`)
- `radius_um`        — radial coordinate at the final time step (µm)
- `radius_normalized`— normalised radius (0 = core, 1 = cell boundary)
- `concentrations`   — `nx × ntr` matrix of final melt-profiles (ppm)

# Element metadata
- `element_names` — `NTuple{N,String}` of element symbols
"""
struct SimulationResult{N}
    time_years::Vector{Float64}
    zircon_radius_um::Vector{Float64}
    growth_rate::Vector{Float64}
    matrix_rate::Vector{Float64}
    rim_conc::Matrix{Float64}
    partition_history::Matrix{Float64}
    interface_melt::Matrix{Float64}
    radius_um::Vector{Float64}
    radius_normalized::Vector{Float64}
    concentrations::Matrix{Float64}
    element_names::NTuple{N,String}
end

"""
    SimulationWorkspace{N}

Pre-allocated workspace for [`simulate!`](@ref).  Create once with
[`SimulationWorkspace(params, elements)`](@ref) and reuse across calls
to avoid all heap allocations.

# Example
```julia
p  = GrowthParams(tfinal_years=500.0)
ed = default_element_data()
ws = SimulationWorkspace(p, ed)
simulate!(ws, 500.0, p, ed)   # zero allocations
```
"""
mutable struct SimulationWorkspace{N}
    # ── grid & temperature ────────────────────────────────────────────────
    xi::Vector{Float64}
    t_phys::Vector{Float64}
    Th::Vector{Float64}
    # ── Zr melt profile + backup ─────────────────────────────────────────
    c::Vector{Float64}
    c0_buf::Vector{Float64}
    # ── trace-element melt profiles ──────────────────────────────────────
    CH0::Matrix{Float64}
    CH_buf::Vector{Float64}
    # ── Thomas-solver workspace ──────────────────────────────────────────
    Av::Vector{Float64}
    Bv::Vector{Float64}
    Cv::Vector{Float64}
    Fv::Vector{Float64}
    sol::Vector{Float64}
    alpha::Vector{Float64}
    beta::Vector{Float64}
    # ── output (written in-place by simulate!) ───────────────────────────
    result::SimulationResult{N}
end

"""
    SimulationWorkspace(params::GrowthParams, elements::ElementData{N}) -> SimulationWorkspace{N}

Allocate all buffers required by [`simulate!`](@ref).
"""
function SimulationWorkspace(params::GrowthParams, elements::ElementData{N}) where {N}
    nt  = params.nt
    nx  = params.nx
    ntr = N
    result = SimulationResult{N}(
        Vector{Float64}(undef, nt),       # time_years
        Vector{Float64}(undef, nt),       # zircon_radius_um
        Vector{Float64}(undef, nt),       # growth_rate
        Vector{Float64}(undef, nt),       # matrix_rate
        Matrix{Float64}(undef, nt, ntr),  # rim_conc
        Matrix{Float64}(undef, nt, ntr),  # partition_history
        Matrix{Float64}(undef, nt, ntr),  # interface_melt
        Vector{Float64}(undef, nx),       # radius_um
        Vector{Float64}(undef, nx),       # radius_normalized
        Matrix{Float64}(undef, nx, ntr),  # concentrations
        elements.names,
    )
    return SimulationWorkspace{N}(
        Vector{Float64}(undef, nx),       # xi
        Vector{Float64}(undef, nt),       # t_phys
        Vector{Float64}(undef, nt),       # Th
        Vector{Float64}(undef, nx),       # c
        Vector{Float64}(undef, nx),       # c0_buf
        Matrix{Float64}(undef, nx, ntr),  # CH0
        Vector{Float64}(undef, nx),       # CH_buf
        Vector{Float64}(undef, nx),       # Av
        Vector{Float64}(undef, nx),       # Bv
        Vector{Float64}(undef, nx),       # Cv
        Vector{Float64}(undef, nx),       # Fv
        Vector{Float64}(undef, nx),       # sol
        Vector{Float64}(undef, nx),       # alpha
        Vector{Float64}(undef, nx),       # beta
        result,
    )
end
