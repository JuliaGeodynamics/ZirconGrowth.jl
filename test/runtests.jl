using Test
using ZirconGrowth

@testset "ZirconGrowth.jl" begin

    # ─── ElementData ──────────────────────────────────────────────────────
    @testset "ElementData loading" begin
        ed = default_element_data()
        @test length(ed) == 11
        @test ed.names[1] == "Hf"
        @test ed.names[2] == "Y"
        @test ed.names[11] == "Li"
        @test ed.concentrations[1] ≈ 6.0
        @test ed.concentrations[8] ≈ 1000.0  # P
    end

    # ─── Zr saturation ───────────────────────────────────────────────────
    @testset "zr_saturation" begin
        # At high T, saturation should be high; at low T, lower
        Csat_high = zr_saturation(1213.0)   # 940 °C
        Csat_low  = zr_saturation(953.0)    # 680 °C
        @test Csat_high > Csat_low
        @test Csat_high > 0.0
        @test Csat_low  > 0.0
        # Known value: at 1213 K with Mfactor=1.3
        # Csat = 490000 / exp(10108/1213 - 1.16*0.3 - 1.48)
        expected = 490_000.0 / exp(10108.0 / 1213.0 - 1.16 * 0.3 - 1.48)
        @test Csat_high ≈ expected rtol=1e-12
    end

    # ─── Zr diffusivity ──────────────────────────────────────────────────
    @testset "zr_diffusivity" begin
        D1 = zr_diffusivity(1213.0, 2.0)
        D2 = zr_diffusivity(953.0, 2.0)
        @test D1 > D2 > 0.0  # diffusion faster at higher T
        # More water → faster diffusion
        D3 = zr_diffusivity(1100.0, 4.0)
        D4 = zr_diffusivity(1100.0, 1.0)
        @test D3 > D4
    end

    # ─── Partition coefficients ──────────────────────────────────────────
    @testset "partition_coefficients" begin
        kd = partition_coefficients(1100.0)
        @test length(kd) == 11
        @test all(k -> k > 0.0, kd)
        # Hf should have highest Kd (it's compatible in zircon)
        @test kd[1] > 100.0  # Hf Kd typically > 100
    end

    # ─── Trace diffusivities ─────────────────────────────────────────────
    @testset "trace_diffusivities" begin
        dif = trace_diffusivities(1100.0, 2.0)
        @test length(dif) == 11
        @test all(d -> d > 0.0, dif)
        # At higher T, diffusion should be faster
        dif_hot  = trace_diffusivities(1200.0, 2.0)
        dif_cold = trace_diffusivities(1000.0, 2.0)
        @test all(dif_hot .> dif_cold)
    end

    # ─── Temperature history ─────────────────────────────────────────────
    @testset "temperature_history" begin
        Td = 1213.0; T0 = 953.0
        # At t=0, T should be Td
        @test temperature_history(Td, T0, 0, 0.0, 0.0, 1.0) ≈ Td
        # At t=tfin, T should be T0
        @test temperature_history(Td, T0, 0, 0.0, 1.0, 1.0) ≈ T0
        # With no oscillations, T is linear and decreasing
        T_mid = temperature_history(Td, T0, 0, 0.0, 0.5, 1.0)
        @test T_mid ≈ 0.5 * (Td + T0)
        # T never drops below T0
        @test temperature_history(Td, T0, 5, 100.0, 0.3, 1.0) >= T0
    end

    # ─── Plagioclase fraction ────────────────────────────────────────────
    @testset "plagioclase_fraction" begin
        @test plagioclase_fraction(1100.0) == 0.0
    end

    # ─── Thomas solver ───────────────────────────────────────────────────
    @testset "thomas_solve!" begin
        # Solve a simple system: diagonal = 2, off-diag = -1
        # i.e. -x_{i-1} + 2x_i - x_{i+1} = F_i
        n = 5
        A = fill(-1.0, n)   # sub-diagonal
        B = fill(-1.0, n)   # super-diagonal
        C = fill(2.0, n)    # diagonal
        F = [1.0, 0.0, 0.0, 0.0, 1.0]
        A[1] = 0.0; B[n] = 0.0  # boundary
        V = zeros(n)
        alpha = zeros(n)
        beta  = zeros(n)
        ZirconGrowth.thomas_solve!(V, A, B, C, F, alpha, beta, n)
        # Check A*x = F
        for i in 1:n
            lhs = C[i] * V[i]
            if i > 1; lhs -= A[i] * V[i-1]; end   # A is sub-diagonal (below)
            if i < n; lhs -= B[i] * V[i+1]; end   # B is super-diagonal (above)
            @test lhs ≈ F[i] atol=1e-12
        end
    end

    # ─── Full simulation (short run) ────────────────────────────────────
    @testset "simulate_zircon_growth — short run" begin
        res = simulate_zircon_growth(500)
        @test length(res.time_years) == res.zircon_radius_um |> length
        @test res.time_years[1] ≈ 0.0
        @test res.time_years[end] ≈ 500.0 rtol=0.01

        # Zircon should grow
        @test res.zircon_radius_um[end] > res.zircon_radius_um[1]

        # Radius should be physically plausible (< 200 µm)
        @test res.zircon_radius_um[end] < 200.0
        @test res.zircon_radius_um[end] > 10.0

        # Rim concentrations should be positive
        @test all(res.rim_conc[end, :] .> 0.0)

        # Check element names are present
        @test length(res.element_names) == 11
        @test "Hf" in res.element_names
        @test "U"  in res.element_names

        # Spatial profiles have right size
        @test size(res.concentrations, 1) == length(res.radius_um)
        @test size(res.concentrations, 2) == 11

        # Normalised radius goes from 0 to ~1
        @test res.radius_normalized[1] ≈ 0.0 atol=1e-10
        @test res.radius_normalized[end] ≈ 1.0 atol=0.01
    end

    # ─── Consistency: compute_profiles == simulate ───────────────────────
    @testset "compute_profiles alias" begin
        r1 = simulate_zircon_growth(500)
        r2 = compute_profiles(500)
        @test r1.zircon_radius_um ≈ r2.zircon_radius_um
        @test r1.rim_conc ≈ r2.rim_conc
    end

    # ─── Longer run produces bigger crystal ──────────────────────────────
    @testset "longer time → bigger crystal" begin
        r_short = simulate_zircon_growth(500)
        p_long = GrowthParams(tfinal_years=1000.0, nt=1000)
        r_long  = simulate_zircon_growth(1000; params=p_long)
        @test r_long.zircon_radius_um[end] > r_short.zircon_radius_um[end]
    end

    # ─── Reproducibility ─────────────────────────────────────────────────
    @testset "deterministic results" begin
        r1 = simulate_zircon_growth(500)
        r2 = simulate_zircon_growth(500)
        @test r1.zircon_radius_um ≈ r2.zircon_radius_um
        @test r1.rim_conc ≈ r2.rim_conc
    end

    # ─── Custom parameters ───────────────────────────────────────────────
    @testset "custom GrowthParams" begin
        p = GrowthParams(tfinal_years=200.0, nx=100, nt=200)
        @test p.nx == 100
        @test p.nt == 200
        res = simulate_zircon_growth(200; params=p)
        @test size(res.concentrations, 1) == 100
        @test length(res.time_years) == 200
    end

    # ─── Edge: very short simulation ─────────────────────────────────────
    @testset "very short simulation" begin
        p = GrowthParams(tfinal_years=10.0, nt=500, nx=200)
        res = simulate_zircon_growth(10; params=p)
        @test res.zircon_radius_um[end] > res.zircon_radius_um[1]
    end

    # ─── Validate against MATLAB reference ───────────────────────────────
    @testset "MATLAB reference comparison" begin
        # The MATLAB single_profile.m with nt=500, nx=200, tfinal=500 yr
        # produces a final zircon radius ≈ 81.455 µm and specific rim
        # concentrations (from growth.csv generated by the MATLAB code).
        p = GrowthParams(tfinal_years=500.0, nt=500, nx=200)
        res = simulate_zircon_growth(500; params=p)
        @test res.zircon_radius_um[end] ≈ 81.4551138255676  rtol=1e-8
        @test res.rim_conc[end, 1]      ≈ 467.12290973884   rtol=1e-6  # Hf
        @test res.rim_conc[end, 2]      ≈ 2936.20715622065  rtol=1e-6  # Y
        @test res.rim_conc[end, 3]      ≈ 89.1206110082121  rtol=1e-6  # U
    end

    # ─── simulate! with pre-allocated workspace ──────────────────────────
    @testset "simulate_zircon_growth! zero-allocation" begin
        ed = default_element_data()
        p  = GrowthParams(tfinal_years=500.0, nt=500, nx=200)
        ws = SimulationWorkspace(p, ed)

        # First call (compiles)
        res = simulate_zircon_growth!(ws, 500.0, p, ed)
        @test res.zircon_radius_um[end] ≈ 81.4551138255676 rtol=1e-8

        # Second call: measure allocations (176 bytes is timer/measurement overhead)
        allocs = @allocated simulate_zircon_growth!(ws, 500.0, p, ed)
        @test allocs <= 176

        # Result matches the allocating wrapper
        res2 = simulate_zircon_growth(500; params=p, elements=ed)
        @test res.zircon_radius_um ≈ res2.zircon_radius_um
        @test res.rim_conc ≈ res2.rim_conc
    end

    # ─── Arbitrary T(t) via Th_override ──────────────────────────────────
    @testset "Th_override: arbitrary non-monotonic T(t)" begin
        # Irregularly-spaced control points in °C and Myr
        # Contains a reheating episode (non-monotonic) to stress-test stability.
        time_Myr = Float64[0.000, 0.005, 0.015, 0.030, 0.040, 0.050, 0.075, 0.150, 0.200]
        T_C      = Float64[940,   920,   870,   810,   840,   860,   790,   720,   680]

        # With custom params (higher resolution)
        p = GrowthParams(time_Myr, T_C; nt=2000, nx=200)
        res = simulate_from_cooling_path(time_Myr, T_C; params=p)

        # All rim concentrations must be finite and positive (no blow-up)
        @test all(isfinite, res.rim_conc)
        @test all(>(0), res.rim_conc)

        # Crystal must have grown from its seed
        @test res.zircon_radius_um[end] > res.zircon_radius_um[1]

        # All time-series outputs finite
        @test all(isfinite, res.zircon_radius_um)
        @test all(isfinite, res.growth_rate)

        # Default params (no explicit GrowthParams needed)
        res2 = simulate_from_cooling_path(time_Myr, T_C)
        @test all(isfinite, res2.rim_conc)
        @test all(>(0), res2.rim_conc)

        # Input validation: mismatched lengths
        @test_throws ArgumentError simulate_from_cooling_path([0.0, 0.1], [900.0])
        # Input validation: unsorted times
        @test_throws ArgumentError simulate_from_cooling_path([0.1, 0.0], [900.0, 800.0])
    end

end
