# ─── Tridiagonal solver and coefficient assembly ──────────────────────────────
# Both functions operate on pre-allocated workspace vectors to avoid allocation.

"""
    thomas_solve!(V, A, B, C, F, alpha, beta, nx)

Solve the tridiagonal system in-place using the Thomas algorithm.

- `A`, `B`, `C`, `F` are the sub-diagonal, super-diagonal, diagonal, and
  right-hand-side vectors (length `nx`).
- `alpha`, `beta` are workspace vectors (length `nx`).
- The solution is written into `V`.

No allocations are performed.
"""
function thomas_solve!(V::AbstractVector{Float64},
                       A::AbstractVector{Float64},
                       B::AbstractVector{Float64},
                       C::AbstractVector{Float64},
                       F::AbstractVector{Float64},
                       alpha::AbstractVector{Float64},
                       beta::AbstractVector{Float64},
                       nx::Int)
    # Forward sweep
    @inbounds alpha[1] = B[1] / C[1]
    @inbounds beta[1]  = F[1] / C[1]
    @inbounds for i in 1:nx-1
        zn = 1.0 / (C[i] - A[i] * alpha[i])
        alpha[i+1] = B[i] * zn
        beta[i+1]  = (F[i] + A[i] * beta[i]) * zn
    end
    # Backward sweep
    @inbounds V[nx] = (F[nx] + A[nx] * beta[nx]) / (C[nx] - A[nx] * alpha[nx])
    @inbounds for i in nx-1:-1:1
        V[i] = V[i+1] * alpha[i+1] + beta[i+1]
    end
    return nothing
end

"""
    assemble_coefficients!(A, B, C, F, c0, xi, H, V, W, R, s, cm,
                           kpl, ktr, nx, Di, dt, itrace)

Fill coefficient vectors `A`, `B`, `C`, `F` for the implicit diffusion
scheme on the non-uniform grid `xi` in the moving-boundary coordinate
system.  `itrace == 0` solves for Zr; `itrace > 0` for the corresponding
trace element.

No allocations are performed; all vectors are modified in-place.
"""
function assemble_coefficients!(A::AbstractVector{Float64},
                                B::AbstractVector{Float64},
                                C_vec::AbstractVector{Float64},
                                F::AbstractVector{Float64},
                                c0::AbstractVector{Float64},
                                xi::AbstractVector{Float64},
                                H::Float64, V::Float64, W::Float64,
                                R::Float64, s::Float64, cm::Float64,
                                kpl::Float64, ktr::Float64,
                                nx::Int, Di::Float64, dt::Float64,
                                itrace::Int)
    # Left boundary (xi = 0, crystal–melt interface)
    if itrace == 0
        # Dirichlet: c = cm (saturation)
        @inbounds A[1] = 0.0
        @inbounds B[1] = 0.0
        @inbounds C_vec[1] = 1.0
        @inbounds F[1] = cm
    else
        # Robin condition for trace elements
        @inbounds C_vec[1] = Di - V * (xi[2] - xi[1]) * (R - s) * (1.0 - ktr)
        @inbounds A[1] = 0.0
        @inbounds B[1] = Di
        @inbounds F[1] = 0.0
    end

    # Interior points
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
        F[i]     = mt * c0[i]
    end

    # Right boundary (xi = 1, cell boundary)
    @inbounds C_vec[nx] = Di + W * (xi[nx] - xi[nx-1]) * (R - s) * (1.0 - kpl)
    @inbounds A[nx] = Di
    @inbounds B[nx] = 0.0
    # F[nx] is not set here — it must be initialised by the caller (typically 0)
    return nothing
end
