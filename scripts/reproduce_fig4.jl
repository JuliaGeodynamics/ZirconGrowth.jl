#!/usr/bin/env julia
"""
    reproduce_fig4.jl

Reproduce Figure 4 from:
  Bindeman & Melnik (2026) Earth Planet. Sci. Lett. 673, 119723.

Figure 4 shows the bulk effective zircon/melt partition coefficient
(Kd_bulk), normalised to the constant input Kd, as a function of log₁₀(Kd)
for elements spanning highly incompatible (Kd = 10⁻³) to highly compatible
(Kd = 10⁴).  Six panels show different crystallisation durations (50, 500,
1000, 5000, 50 000, 500 000 yr).  Within each panel, three curves represent
trace elements diffusing slower (0.25×), equal (1×), or faster (4×) than Zr.

Default parameters: 940→680 °C, 2 wt% H₂O, M = 1.3.

Run from the package root:
    julia --project scripts/reproduce_fig4.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ZirconGrowth
using GLMakie

# ─── Lower-level imports from ZirconGrowth ────────────────────────────────────
# We need the solver internals for the parametric trace-element sweep.
using ZirconGrowth: assemble_coefficients!, thomas_solve!,
                    zr_saturation, zr_diffusivity, temperature_history,
                    plagioclase_fraction_poly, SECONDS_PER_YEAR

"""
    bulk_kd_sweep(tfinal_years, kd_values, d_factors; kwargs...)

For a given crystallisation duration, compute the bulk effective Kd for each
combination of constant partition coefficient (`kd_values`) and diffusion
scaling factor (`d_factors`, relative to D_Zr).

Returns a matrix of size `(length(kd_values), length(d_factors))` containing
Kd_bulk / Kd for each combination.
"""
function bulk_kd_sweep(tfinal_years::Float64,
                       kd_values::Vector{Float64},
                       d_factors::Vector{Float64};
                       nt::Int = clamp(round(Int, tfinal_years), 500, 5000),
                       nx::Int = 300,
                       Td::Float64 = 940.0 + 273.0,
                       T0::Float64 = 680.0 + 273.0,
                       Rc::Float64 = 0.2,
                       s0::Float64 = 1e-4,
                       Czirc::Float64 = 490_000.0,
                       kpl::Float64 = 0.2,
                       XH2O::Float64 = 2.0,
                       Mfac::Float64 = 1.3)

    tyear = SECONDS_PER_YEAR
    tfin  = tfinal_years * tyear

    # ── Grid ─────────────────────────────────────────────────────────────
    xi = [1.0 - cos(Float64(pi) * (i - 1) / (nx - 1) / 2.0)^2 for i in 1:nx]

    # ── Time & temperature ───────────────────────────────────────────────
    dt_phys = tfin / (nt - 1)
    t_phys  = [(i - 1) * dt_phys for i in 1:nt]
    Th      = [temperature_history(Td, T0, 0, 0.0, t_phys[i], tfin) for i in 1:nt]

    # ── Scaling ──────────────────────────────────────────────────────────
    Dscale = zr_diffusivity(0.5 * (Th[1] + Th[nt]), 3.0)
    tscale = Rc^2 / Dscale
    dt     = dt_phys / tscale
    s      = s0 / Rc
    R      = 1.0

    # ── Phase 1: Solve Zr growth (store per-step geometry) ───────────────
    c = fill(zr_saturation(Td; Mfactor=Mfac, Czirc=Czirc), nx)

    # Storage for geometry at each step
    s_hist = Vector{Float64}(undef, nt)
    R_hist = Vector{Float64}(undef, nt)
    V_hist = Vector{Float64}(undef, nt)
    W_hist = Vector{Float64}(undef, nt)
    H_hist = Vector{Float64}(undef, nt)
    cm_hist = Vector{Float64}(undef, nt)
    DZr_hist = Vector{Float64}(undef, nt)

    s_hist[1] = s;  R_hist[1] = R;  V_hist[1] = 0.0;  W_hist[1] = 0.0
    H_hist[1] = R - s;  cm_hist[1] = c[1]
    DZr_hist[1] = zr_diffusivity(Th[1], XH2O) / Dscale

    # Workspace vectors for Thomas solver (Zr)
    Av = zeros(nx);  Bv = zeros(nx);  Cv = zeros(nx);  Fv = zeros(nx)
    sol = zeros(nx);  alpha = zeros(nx);  beta = zeros(nx)
    c0_buf = zeros(nx)

    for j in 2:nt
        copyto!(c0_buf, c)
        s0_loc = s
        R0     = R
        errC   = 1e10
        itt    = 0

        Di_zr = zr_diffusivity(Th[j], XH2O) / Dscale
        cm    = zr_saturation(Th[j]; Mfactor=Mfac, Czirc=Czirc)
        dXpdT = (plagioclase_fraction_poly(Th[j] + 1e-5) -
                 plagioclase_fraction_poly(Th[j] - 1e-5)) / 2e-5
        dTdt  = (Th[j] - Th[j-1]) / dt

        V_loc = 0.0;  W_loc = 0.0;  H_loc = 0.0
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

            errC = maximum(abs(sol[k] - c[k]) for k in 1:nx)
            copyto!(c, sol)
            itt += 1
        end

        s_hist[j] = s;  R_hist[j] = R;  V_hist[j] = V_loc;  W_hist[j] = W_loc
        H_hist[j] = H_loc;  cm_hist[j] = cm;  DZr_hist[j] = Di_zr
    end

    # ── Phase 2: For each (Kd, D_factor), solve trace element & compute bulk Kd
    nkd = length(kd_values)
    ndf = length(d_factors)
    result = Matrix{Float64}(undef, nkd, ndf)

    # workspace for trace solve
    CH0   = zeros(nx)
    CH_buf = zeros(nx)
    At = zeros(nx);  Bt = zeros(nx);  Ct = zeros(nx);  Ft = zeros(nx)
    sol_t = zeros(nx);  alpha_t = zeros(nx);  beta_t = zeros(nx)

    C0_melt = 100.0  # arbitrary initial melt concentration (ppm)

    for (ikd, kd_const) in enumerate(kd_values)
        for (idf, dfac) in enumerate(d_factors)
            # Reset trace element melt profile to uniform
            fill!(CH0, C0_melt)

            # Track zircon concentration history for bulk averaging
            # Bulk Kd = (mass-weighted avg C in zircon) / (volume-weighted avg C in melt)
            # Zircon grows shell by shell: each timestep adds a shell of concentration kd_const * CH0[1]
            zircon_mass_sum = 0.0  # ∑ C_zircon_i * ΔV_i
            zircon_vol_sum  = 0.0  # ∑ ΔV_i

            # Initial shell
            s_prev = s_hist[1]
            zircon_mass_sum += kd_const * C0_melt * (4π / 3) * s_prev^3
            zircon_vol_sum  += (4π / 3) * s_prev^3

            for j in 2:nt
                ktr   = kd_const
                Di_el = DZr_hist[j] * dfac  # scale D relative to Zr

                Ft[nx] = 0.0
                assemble_coefficients!(At, Bt, Ct, Ft, CH0, xi,
                                       H_hist[j], V_hist[j], W_hist[j],
                                       R_hist[j], s_hist[j], cm_hist[j],
                                       0.0, ktr, nx, Di_el, dt, 1)
                thomas_solve!(sol_t, At, Bt, Ct, Ft, alpha_t, beta_t, nx)
                copyto!(CH0, sol_t)

                # Shell volume: 4π/3 * (s_j³ - s_{j-1}³)
                s_cur  = s_hist[j]
                dV = (4π / 3) * abs(s_cur^3 - s_prev^3)
                c_zircon_shell = ktr * CH0[1]  # instantaneous interface concentration × Kd
                zircon_mass_sum += c_zircon_shell * dV
                zircon_vol_sum  += dV
                s_prev = s_cur
            end

            # Bulk melt: volume-weighted average over final melt profile
            # Use trapezoidal integration on the physical radius grid
            r_phys = [(s_hist[nt] + H_hist[nt] * xi[i]) * Rc for i in 1:nx]
            melt_mass = 0.0
            melt_vol  = 0.0
            for i in 2:nx
                dr = r_phys[i] - r_phys[i-1]
                r_mid = 0.5 * (r_phys[i] + r_phys[i-1])
                c_mid = 0.5 * (CH0[i] + CH0[i-1])
                shell_vol = 4π * r_mid^2 * dr
                melt_mass += c_mid * shell_vol
                melt_vol  += shell_vol
            end

            C_zircon_bulk = zircon_vol_sum > 0.0 ? zircon_mass_sum / zircon_vol_sum : 0.0
            C_melt_bulk   = melt_vol > 0.0 ? melt_mass / melt_vol : C0_melt

            Kd_bulk = C_melt_bulk > 0.0 ? C_zircon_bulk / C_melt_bulk : 0.0
            result[ikd, idf] = Kd_bulk / kd_const  # normalised
        end
    end

    return result
end

function run_fig4()
    # ── Parametric sweep parameters ──────────────────────────────────────
    log_kd_range = range(-3, 4, length=50)
    kd_values    = 10.0 .^ collect(log_kd_range)

    d_factors = [0.25, 1.0, 4.0]
    d_labels  = ["D = 0.25 × D_Zr", "D = D_Zr", "D = 4 × D_Zr"]
    d_colors  = [:blue, :black, :red]

    durations = [50.0, 500.0, 1000.0, 5000.0, 50_000.0, 500_000.0]
    dur_labels = ["50 yr", "500 yr", "1000 yr", "5000 yr", "50 000 yr", "500 000 yr"]

    fig = Figure(size=(1400, 900))

    for (ipanel, (tfyr, dur_label)) in enumerate(zip(durations, dur_labels))
        row = (ipanel - 1) ÷ 3 + 1
        col = (ipanel - 1) % 3 + 1

        @info "Running sweep for $dur_label…"
        result = bulk_kd_sweep(tfyr, kd_values, d_factors)

        ax = Axis(fig[row, col];
                  xlabel = row == 2 ? "log₁₀(Kd)" : "",
                  ylabel = col == 1 ? "Kd_bulk / Kd" : "",
                  title  = dur_label,
                  yscale = log10,
                  limits = ((-3, 4), (1e-2, 1e4)))

        hlines!(ax, [1.0]; color=:gray, linestyle=:dash, linewidth=0.5)

        for (idf, (lbl, clr)) in enumerate(zip(d_labels, d_colors))
            y = result[:, idf]
            # Clamp for log scale display
            y = max.(y, 1e-3)
            lines!(ax, collect(log_kd_range), y;
                   linewidth=2, color=clr, label=lbl)
        end

        if ipanel == 1
            axislegend(ax; position=:lt, framevisible=false)
        end
    end

    return fig
end

fig = run_fig4()
display(fig)
@info "Figure 4 displayed."
