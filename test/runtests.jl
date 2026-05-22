using Test
using PlasticViz

@testset "PlasticViz Mohr-Coulomb helpers" begin

    @testset "mohr_failure_state — safe" begin
        # c=10, φ=30°, σ₁=20, σ₃=5
        # P=12.5, R=7.5, tanφ≈0.577, d=(10+12.5*0.577)/1.155≈14.9 > R → safe
        state, d, tp, θ = PlasticViz.mohr_failure_state(10.0, 30.0, 20.0, 5.0)
        @test state == :safe
        @test d > 7.5
        @test θ ≈ 60.0
        @test tp[2] > 0.0  # tangent point is in upper half-plane
    end

    @testset "mohr_failure_state — critical (well past envelope)" begin
        # c=5, φ=20°, σ₁=100, σ₃=0
        # P=50, R=50, d≈21.8 — circle far outside envelope → :critical
        state, d, _, _ = PlasticViz.mohr_failure_state(5.0, 20.0, 100.0, 0.0)
        @test state == :critical
        @test d < 50.0
    end

    @testset "mohr_failure_state — failure angle" begin
        _, _, _, θ = PlasticViz.mohr_failure_state(15.0, 25.0, 50.0, 10.0)
        @test θ ≈ 57.5
    end

    @testset "mohr_failure_state — tangent point on circle" begin
        # Tangent point must lie exactly on the Mohr circle
        c, φ, σ₁, σ₃ = 15.0, 30.0, 80.0, 10.0
        _, _, tp, _ = PlasticViz.mohr_failure_state(c, φ, σ₁, σ₃)
        P = (σ₁ + σ₃) / 2.0
        R = (σ₁ - σ₃) / 2.0
        dist_from_centre = sqrt((tp[1] - P)^2 + tp[2]^2)
        @test dist_from_centre ≈ R atol=1e-4
    end

    @testset "failure_block_geometry — failure angle" begin
        geo = PlasticViz.failure_block_geometry(30.0, 80.0, 10.0)
        @test geo.θ ≈ 60.0
    end

    @testset "failure_block_geometry — plane endpoints on rectangle boundary" begin
        geo = PlasticViz.failure_block_geometry(30.0, 80.0, 10.0)
        p1, p2 = geo.plane_pts
        on_boundary = p -> isapprox(abs(p[1]), 1.0, atol=1e-5) ||
                           isapprox(abs(p[2]), 1.0, atol=1e-5)
        @test on_boundary(p1)
        @test on_boundary(p2)
    end

    @testset "failure_block_geometry — arrow count" begin
        geo = PlasticViz.failure_block_geometry(30.0, 80.0, 10.0)
        @test length(geo.σ₁_tails) == 6  # 3 top + 3 bottom
        @test length(geo.σ₁_vecs)  == 6
    end

    @testset "failure_block_geometry — zero σ₃ hides σ₃ arrows" begin
        geo = PlasticViz.failure_block_geometry(30.0, 80.0, 0.0)
        @test !geo.show_σ₃
    end

    @testset "mohr_failure_state — critical" begin
        # c=10, φ=0°: d=c=10, R=(σ₁-σ₃)/2. Set R=10 → σ₁=30, σ₃=10
        state, _, _, _ = PlasticViz.mohr_failure_state(10.0, 0.0, 30.0, 10.0)
        @test state == :critical
    end

    @testset "mohr_failure_state — tensile failure" begin
        # T₀=5 MPa, σ₃=-10 → σ₃ ≤ -T₀+tol → tensile failure regardless of shear envelope
        state, _, _, _ = PlasticViz.mohr_failure_state(20.0, 30.0, 50.0, -10.0; T₀=5.0)
        @test state == :tensile
    end

    @testset "mohr_failure_state — tensile at boundary" begin
        # σ₃ = -T₀ exactly → tensile (slider sits at the floor, ≤ check with tol triggers)
        state, _, _, _ = PlasticViz.mohr_failure_state(10.0, 30.0, 50.0, -5.0; T₀=5.0)
        @test state == :tensile
    end

    @testset "mohr_failure_state — T₀=Inf never triggers tensile" begin
        # Default T₀=Inf: even very negative σ₃ should not give :tensile
        state, _, _, _ = PlasticViz.mohr_failure_state(20.0, 30.0, 100.0, -50.0)
        @test state != :tensile
    end

end
