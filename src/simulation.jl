# ─── Main simulation loop ─────────────────────────────────────────────────────

"""
    simulate_zircon_growth!(ws::SimulationWorkspace{N}, tfinal_years, params, elements) -> SimulationResult{N}

Run the zircon growth model **in-place** into the pre-allocated workspace `ws`.
Performs **zero heap allocations** — all temporary and output arrays live inside
`ws`.  Returns `ws.result` for convenience.

Create the workspace once with [`SimulationWorkspace(params, elements)`](@ref):

```julia
ws = SimulationWorkspace(params, elements)
simulate_zircon_growth!(ws, 500.0, params, elements)   # 0 allocations
simulate_zircon_growth!(ws, 500.0, params, elements)   # reuse — still 0 allocations
```
"""
function simulate_zircon_growth!(ws::SimulationWorkspace{N},
                   tfinal_years::Float64,
                   params::GrowthParams,
                   elements::ElementData{N}) where {N}

    nt    = params.nt
    nx    = params.nx
    Rc    = params.Rc
    s0    = params.s0
    Td    = params.Td
    T0    = params.T0
    Czirc = params.Czirc
    kpl   = params.kpl
    XH2O  = params.XH2O
    Mfac  = params.Mfactor
    n_osc = params.n_oscillations
    amp   = params.amplitude
    use_mm = params.major_minerals

    tyear = SECONDS_PER_YEAR
    tfin  = tfinal_years * tyear

    # Aliases into workspace
    xi     = ws.xi
    t_phys = ws.t_phys
    Th     = ws.Th
    c      = ws.c
    c0_buf = ws.c0_buf
    CH0    = ws.CH0
    CH_buf = ws.CH_buf
    Av     = ws.Av
    Bv     = ws.Bv
    Cv     = ws.Cv
    Fv     = ws.Fv
    sol    = ws.sol
    alpha  = ws.alpha
    beta   = ws.beta

    res         = ws.result
    Xs          = res.zircon_radius_um
    V_hist      = res.growth_rate
    W_hist      = res.matrix_rate
    CHS         = res.rim_conc
    partHist    = res.partition_history
    ifaceMelt   = res.interface_melt
    time_years  = res.time_years
    radius_um   = res.radius_um
    radius_norm = res.radius_normalized
    concOut     = res.concentrations

    ntr = N

    # Cosine-clustered grid
    @inbounds for i in 1:nx
        xi[i] = 1.0 - cos(Float64(pi) * (i - 1) / (nx - 1) / 2.0)^2
    end

    # Time array (physical seconds)
    dt_phys = tfin / (nt - 1)
    @inbounds for i in 1:nt
        t_phys[i] = (i - 1) * dt_phys
    end

    # Temperature at every step
    @inbounds for i in 1:nt
        Th[i] = temperature_history(Td, T0, n_osc, amp, t_phys[i], tfin)
    end

    # Trace-element melt profiles (uniform initial)
    @inbounds for j in 1:ntr
        c0j = elements.concentrations[j]
        for i in 1:nx
            CH0[i, j] = c0j
        end
    end

    # Zero output matrices
    fill!(CHS, 0.0)
    fill!(partHist, 0.0)
    fill!(ifaceMelt, 0.0)
    fill!(Xs, 0.0)
    fill!(V_hist, 0.0)
    fill!(W_hist, 0.0)

    kt0 = partition_coefficients(Th[1]; Mfactor=Mfac, Czirc=Czirc)
    @inbounds for j in 1:ntr
        CHS[1, j]      = elements.concentrations[j] * kt0[j]
        partHist[1, j]  = kt0[j]
        ifaceMelt[1, j] = elements.concentrations[j]
    end

    # Zr melt profile (uniform at saturation)
    csat0 = zr_saturation(Td; Mfactor=Mfac, Czirc=Czirc)
    @inbounds for i in 1:nx
        c[i] = csat0
    end

    # Dimensionless scaling
    Dscale = zr_diffusivity(0.5 * (Th[1] + Th[nt]), 3.0)
    tscale = Rc^2 / Dscale
    dt     = dt_phys / tscale
    s      = s0 / Rc
    R      = 1.0
    Xs[1]  = s * Rc * 1e4
    V0     = Rc * 1e4 / tscale * tyear

    # Main time loop (zero allocations)
    Mpl = 0.0
    @inbounds for j in 2:nt
        for k in 1:nx; c0_buf[k] = c[k]; end
        s0_loc = s
        R0     = R
        errC   = 1e10
        itt    = 0

        Di_zr = zr_diffusivity(Th[j], XH2O) / Dscale
        cm    = zr_saturation(Th[j]; Mfactor=Mfac, Czirc=Czirc)

        dXpdT = if use_mm
            (plagioclase_fraction_poly(Th[j] + 1e-5) -
             plagioclase_fraction_poly(Th[j] - 1e-5)) / 2e-5
        else
            0.0
        end
        dTdt  = (Th[j] - Th[j-1]) / dt

        # Picard iteration for Zr growth
        V_loc = 0.0
        W_loc = 0.0
        H_loc = 0.0
        while errC > 1e-7 && itt < 1000
            V_loc = -Di_zr * (c[2] - cm) / (xi[2] - xi[1]) / (R - s) / (cm - Czirc)
            W_loc = -dXpdT * dTdt / 3.0 / R^2
            R = R0 + W_loc * dt
            s = s0_loc + V_loc * dt
            H_loc = R - s

            Fv[nx] = 0.0
            assemble_coefficients!(Av, Bv, Cv, Fv, c0_buf, xi,
                                   H_loc, V_loc, W_loc, R, s, cm,
                                   kpl, 0.0, nx, Di_zr, dt, 0)
            thomas_solve!(sol, Av, Bv, Cv, Fv, alpha, beta, nx)

            errC = 0.0
            for k in 1:nx
                d = abs(sol[k] - c[k])
                errC = d > errC ? d : errC
            end
            for k in 1:nx; c[k] = sol[k]; end
            itt += 1
        end

        # Trace elements (one linear solve each)
        kt  = partition_coefficients(Th[j]; Mfactor=Mfac, Czirc=Czirc)
        Dif = trace_diffusivities(Th[j], XH2O)
        for el in 1:ntr
            partHist[j, el] = kt[el]
        end

        for el in 1:ntr
            ktr   = kt[el]
            Di_el = Dif[el] / Dscale

            Fv[nx] = 0.0
            _assemble_col!(Av, Bv, Cv, Fv, CH0, el, xi,
                           H_loc, V_loc, W_loc, R, s, cm,
                           0.0, ktr, nx, Di_el, dt, el)
            thomas_solve!(CH_buf, Av, Bv, Cv, Fv, alpha, beta, nx)

            ifaceMelt[j, el] = CH_buf[1]
            CHS[j, el]       = ktr * CH_buf[1]
            for k in 1:nx
                CH0[k, el] = CH_buf[k]
            end
        end

        Mpl += -4.0 * Float64(pi) * R^2 * W_loc * dt * kpl * c[nx]
        Xs[j]     = s * Rc * 1e4
        V_hist[j] = V_loc
        W_hist[j] = W_loc
    end

    # Build output arrays (in-place)
    @inbounds for i in 1:nt
        time_years[i] = t_phys[i] / tyear
        V_hist[i]     = V_hist[i] * V0
        W_hist[i]     = W_hist[i] * V0
    end

    @inbounds for i in 1:nx
        radius_um[i] = (s + (R - s) * xi[i]) * Rc * 1e4
    end
    rspan = radius_um[nx] - radius_um[1]
    if abs(rspan) < eps(_absmax(radius_um, nx) + 1.0)
        @inbounds for i in 1:nx; radius_norm[i] = 0.0; end
    else
        inv_rspan = 1.0 / rspan
        r1 = radius_um[1]
        @inbounds for i in 1:nx
            radius_norm[i] = (radius_um[i] - r1) * inv_rspan
        end
    end

    # Copy final melt profiles into result
    @inbounds for j in 1:ntr
        for i in 1:nx
            concOut[i, j] = CH0[i, j]
        end
    end

    return res
end

@inline function _absmax(v::AbstractVector{Float64}, n::Int)
    m = 0.0
    @inbounds for i in 1:n
        a = abs(v[i])
        m = a > m ? a : m
    end
    return m
end

# Assemble coefficients reading directly from a matrix column (avoids view allocation)
function _assemble_col!(A::Vector{Float64}, B::Vector{Float64},
                        C_vec::Vector{Float64}, F::Vector{Float64},
                        CH0::Matrix{Float64}, col::Int,
                        xi::Vector{Float64},
                        H::Float64, V::Float64, W::Float64,
                        R::Float64, s::Float64, cm::Float64,
                        kpl::Float64, ktr::Float64,
                        nx::Int, Di::Float64, dt::Float64,
                        itrace::Int)
    @inbounds C_vec[1] = Di - V * (xi[2] - xi[1]) * (R - s) * (1.0 - ktr)
    @inbounds A[1] = 0.0
    @inbounds B[1] = Di
    @inbounds F[1] = 0.0

    @inbounds for i in 2:nx-1
        xi12 = 0.5 * (xi[i-1] + xi[i])
        h12  = xi[i] - xi[i-1]
        xi23 = 0.5 * (xi[i+1] + xi[i])
        h23  = xi[i+1] - xi[i]

        t3  = H * (V - W)
        t4  = xi12 * xi12
        t10 = (2.0 * V - W) * s - V * R
        t12 = V * s
        t25 = xi23 * xi23

        A[i] = -(dt / h12 * (h12 * H * (xi12 * t10 + t4 * t3 - Di - t12) +
                 2.0 * (xi12 * H + s) * Di)) * 0.5
        B[i] =  (dt / h23 * (h23 * H * (xi23 * t10 + t25 * t3 - Di - t12) -
                 2.0 * (xi23 * H + s) * Di)) * 0.5

        H2 = H * H
        mt = 0.5 * ((-xi12 + 2.0 - xi23) * s + R * (xi12 + xi23)) * H2 * (xi12 - xi23)

        C_vec[i] = A[i] + B[i] + mt
        F[i]     = mt * CH0[i, col]
    end

    @inbounds C_vec[nx] = Di + W * (xi[nx] - xi[nx-1]) * (R - s) * (1.0 - kpl)
    @inbounds A[nx] = Di
    @inbounds B[nx] = 0.0
    return nothing
end

# ─── Convenience wrapper (allocating) ─────────────────────────────────────────

"""
    simulate_zircon_growth(tfinal_years; params, elements) -> SimulationResult

Allocating convenience wrapper around [`simulate_zircon_growth!`](@ref).  Allocates a
fresh [`SimulationWorkspace`](@ref), runs the simulation, and returns the
result.  For zero-allocation repeated runs, use `simulate_zircon_growth!` with a
pre-allocated workspace instead.
"""
function simulate_zircon_growth(tfinal_years::Real;
                  params::GrowthParams = GrowthParams(tfinal_years=Float64(tfinal_years)),
                  elements::ElementData{N} = default_element_data()) where {N}
    ws = SimulationWorkspace(params, elements)
    return simulate_zircon_growth!(ws, Float64(tfinal_years), params, elements)
end

"""
    compute_profiles(tfinal_years; kwargs...) -> SimulationResult

Alias for [`simulate_zircon_growth`](@ref).
"""
compute_profiles(tfinal_years::Real; kwargs...) = simulate_zircon_growth(tfinal_years; kwargs...)
