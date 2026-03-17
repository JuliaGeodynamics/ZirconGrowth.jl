# ─── Temperature-dependent physical properties ────────────────────────────────
# All functions are written to be allocation-free (scalar in, scalar out).

const SECONDS_PER_YEAR = 365.0 * 3600.0 * 24.0

# ── Zr saturation in the melt ─────────────────────────────────────────────────

"""
    zr_saturation(T; Mfactor=1.3, Czirc=490_000.0) -> Float64

Zr saturation concentration (ppm) in the melt at temperature `T` (K),
following Boehnke et al. (2013) Chem. Geol. 351, 324.

```math
C_{\\text{sat}} = \\frac{C_{\\text{zirc}}}{\\exp\\!\\bigl(10108/T - 1.16(M-1) - 1.48\\bigr)}
```
"""
@inline function zr_saturation(T::Float64; Mfactor::Float64=1.3, Czirc::Float64=490_000.0)
    return Czirc / exp(10108.0 / T - 1.16 * (Mfactor - 1.0) - 1.48)
end

# ── Zr diffusion in the melt ──────────────────────────────────────────────────

"""
    zr_diffusivity(T, XH2O) -> Float64

Zr diffusion coefficient (cm² s⁻¹) in silicate melt at temperature `T` (K)
with `XH2O` wt% dissolved water.  Parameterisation from the best-fit with
H₂O dependence.

```math
\\ln D = -\\frac{11.4x + 3.13}{0.84x + 1}
        -\\frac{21.4x + 47}{1.06x + 1}\\,\\theta,
\\qquad \\theta = 1000/T
```

Returns `D` in cm² s⁻¹ (the factor 10⁴ converts from m² s⁻¹ to cm² s⁻¹).
"""
@inline function zr_diffusivity(T::Float64, XH2O::Float64)
    θ = 1000.0 / T
    lnD = -(11.4 * XH2O + 3.13) / (0.84 * XH2O + 1.0) -
          (21.4 * XH2O + 47.0) / (1.06 * XH2O + 1.0) * θ
    return exp(lnD) * 1e4   # cm²/s
end

# ── Trace-element partition coefficients ───────────────────────────────────────

"""
    partition_coefficients(T; Mfactor=1.3, Czirc=490_000.0) -> NTuple{11,Float64}

Temperature-dependent zircon/melt partition coefficients ``K_d`` for the 11
trace elements (Hf, Y, U, Th, Sm, Dy, Yb, P, Nb, Sc, Li) at temperature
`T` (K).  Each ``\\ln K_d`` is a linear or saturation-dependent function of
``\\theta = 1000/T``.
"""
@inline function partition_coefficients(T::Float64;
                                        Mfactor::Float64=1.3,
                                        Czirc::Float64=490_000.0)
    X = 1000.0 / T
    Csat = zr_saturation(T; Mfactor, Czirc)

    lnKd_Hf = 11.29 * X - 2.275
    lnKd_Y  = 19.47 * X - 13.04
    lnKd_U  = 15.32 * X - 9.17
    lnKd_Th = 13.02 * X - 8.354
    lnKd_Sm = log(13.338 * Csat^(-0.622))
    lnKd_Dy = log(2460.0 * Csat^(-0.867))
    lnKd_Yb = log(33460.0 * Csat^(-1.040))
    lnKd_P  =  7.646 * X - 5.047
    lnKd_Nb = 11.29 * X - 11.626
    lnKd_Sc = 19.47 * X - 19.6
    lnKd_Li = 11.29 * X - 13.625

    return (exp(lnKd_Hf), exp(lnKd_Y), exp(lnKd_U), exp(lnKd_Th),
            exp(lnKd_Sm), exp(lnKd_Dy), exp(lnKd_Yb), exp(lnKd_P),
            exp(lnKd_Nb), exp(lnKd_Sc), exp(lnKd_Li))
end

# ── Trace-element diffusion coefficients ───────────────────────────────────────

"""
    trace_diffusivities(T, XH2O) -> NTuple{11,Float64}

Diffusion coefficients (cm² s⁻¹) for the 11 trace elements in silicate melt
at temperature `T` (K) and water content `XH2O` (wt%).  Pressure is fixed
at 0.1 GPa.
"""
@inline function trace_diffusivities(T::Float64, XH2O::Float64)
    P  = 0.1  # GPa
    θ  = 1000.0 / T
    x  = XH2O
    sx = sqrt(x)
    log1e4 = log(1e4)

    # Hf: same as Zr / 1e5 (note: MATLAB divides D_Zr by 10e4 = 1e5)
    lnD_Hf = log(zr_diffusivity(T, x) / 1e5)

    lnD_Y  = -9.925340370 - (35587.24428 - 2615.214492 * x) / T
    lnD_U  = -6.37 - 2.65 * sx - (44729.0 - 1093.0 * P - 8944.0 * sx) / T
    lnD_Th = -7.02 - 1.30 * sx - (44682.0 - 1370.0 * P - 7281.0 * sx) / T
    lnD_P  = (9.469 * x - 181.5) / (x + 2.374) * (θ - 0.54) - 26.7
    lnD_Nb = -15.77 - 21329.0 / T - log1e4
    lnD_Sc = -12.88 - 22717.0 / T - log1e4
    lnD_Li = -13.04 - 10966.0 / T - log1e4
    lnD_Yb = -15.68 - 22478.0 / T - log1e4
    lnD_Dy = -1.93  - 43180.0 / T - log1e4
    lnD_Sm = -5.91  - 39377.0 / T - log1e4

    f = 1e4  # m²/s → cm²/s
    return (exp(lnD_Hf) * f, exp(lnD_Y) * f, exp(lnD_U) * f, exp(lnD_Th) * f,
            exp(lnD_Sm) * f, exp(lnD_Dy) * f, exp(lnD_Yb) * f, exp(lnD_P) * f,
            exp(lnD_Nb) * f, exp(lnD_Sc) * f, exp(lnD_Li) * f)
end

# ── Temperature history ────────────────────────────────────────────────────────

"""
    temperature_history(Td, T0, n, amp, t) -> Float64

Temperature (K) at time `t` given a linear cooling path from `Td` to `T0`
with optional sinusoidal oscillations (`n` cycles, amplitude `amp` in K).
"""
@inline function temperature_history(Td::Float64, T0::Float64,
                                     n::Int, amp::Float64,
                                     t::Float64, tfin::Float64)
    return max(T0, Td - (Td - T0) * t / tfin - amp * sin(n * 2.0 * Float64(pi) * t / tfin))
end

# ── Plagioclase fraction ──────────────────────────────────────────────────────

"""
    plagioclase_fraction(T) -> Float64

Volume fraction of plagioclase crystallised at temperature `T` (K).
Returns 0 (disabled in the reference MATLAB implementation).
"""
@inline function plagioclase_fraction(::Float64)
    return 0.0
end

"""
    plagioclase_fraction_poly(T) -> Float64

Volume fraction of plagioclase crystallised at temperature `T` (K),
using the 4th-order polynomial fit from the MATLAB implementation.
Used when `major_minerals = true` in [`GrowthParams`](@ref).
"""
@inline function plagioclase_fraction_poly(T::Float64)
    x = T - 273.0
    Xpl = 1.94319815919224e-07 * x^4 - 6.44583169395759e-04 * x^3 +
          0.799105630147523 * x^2 - 439.159560503417 * x + 90385.6027471822
    return Xpl / 100.0
end
