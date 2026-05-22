# Helper: classify a Mohr-Coulomb stress state and return diagnostic geometry.
# Returns: (state::Symbol, d::Float64, tangent_pt::Point2f, θ::Float64)
#   state       — :safe | :critical | :failed
#   d           — perpendicular distance from circle centre to failure line
#   tangent_pt  — closest point on circle to the failure line
#   θ           — failure plane angle in degrees (45 + φ/2)
function mohr_failure_state(c::Real, φ::Real, σ₁::Real, σ₃::Real)
    P    = (σ₁ + σ₃) / 2.0
    R    = (σ₁ - σ₃) / 2.0
    tanφ = tand(φ)
    norm = sqrt(1.0 + tanφ^2)
    d    = (c + P * tanφ) / norm          # perp. distance centre → line
    tol  = 1e-2
    state = R > d + tol ? :failed : (R > d - tol ? :critical : :safe)
    θ     = 45.0 + φ / 2.0
    # Tangent point = centre + R * unit_normal_toward_failure_line
    # Line: -tanφ·σ + τ - c = 0  →  gradient (-tanφ, 1), pointing from sub- to super-yield
    tangent_pt = Point2f(P - R * tanφ / norm, R / norm)
    return state, d, tangent_pt, θ
end

# Helper: compute all drawable geometry for the block sketch panel.
# All coordinates live in normalised [-1,1]×[-1,1] block space.
# Returns a NamedTuple with:
#   plane_pts   — Tuple{Point2f, Point2f} clipped to rectangle boundary
#   θ           — failure plane angle (degrees)
#   σ₁_tails, σ₁_vecs — arrows for σ₁ (vertical, compressive)
#   σ₃_tails, σ₃_vecs — arrows for σ₃ (horizontal, compressive); empty when σ₃ ≈ 0
function failure_block_geometry(φ::Real, σ₁::Real, σ₃::Real)
    θ    = 45.0 + φ / 2.0
    tanθ = tand(θ)

    # Clip line y = tanθ·x to rectangle [-1,1]×[-1,1]
    x_at_ytop = abs(tanθ) > 1e-9 ? 1.0f0 / Float32(tanθ) : Inf32
    p1, p2 = if isfinite(x_at_ytop) && abs(x_at_ytop) ≤ 1.0f0
        Point2f(-x_at_ytop, -1f0), Point2f(x_at_ytop, 1f0)
    else
        y_r = clamp(Float32(tanθ),  -1f0, 1f0)
        y_l = clamp(Float32(-tanθ), -1f0, 1f0)
        Point2f(-1f0, y_l), Point2f(1f0, y_r)
    end

    # Arrow lengths scale with stress magnitude (σ₁ always ≥ minimum visible)
    σ₁_len = clamp(0.15f0 + 0.4f0 * Float32(σ₁) / 150f0, 0.15f0, 0.55f0)
    σ₃_len = 0.15f0 + 0.4f0 * Float32(σ₃) / 150f0

    xs = (-0.5f0, 0f0, 0.5f0)
    σ₁_tails = vcat([Point2f(x,  1f0 + σ₁_len) for x in xs],
                    [Point2f(x, -1f0 - σ₁_len) for x in xs])
    σ₁_vecs  = vcat([Vec2f(0f0, -σ₁_len) for _ in xs],
                    [Vec2f(0f0,  σ₁_len) for _ in xs])

    if σ₃ < 0.5
        σ₃_tails = Point2f[]
        σ₃_vecs  = Vec2f[]
    else
        ys = (-0.4f0, 0f0, 0.4f0)
        σ₃_tails = vcat([Point2f( 1f0 + σ₃_len, y) for y in ys],
                        [Point2f(-1f0 - σ₃_len, y) for y in ys])
        σ₃_vecs  = vcat([Vec2f(-σ₃_len, 0f0) for _ in ys],
                        [Vec2f( σ₃_len, 0f0) for _ in ys])
    end

    return (plane_pts=(p1, p2), θ=θ,
            σ₁_tails=σ₁_tails, σ₁_vecs=σ₁_vecs,
            σ₃_tails=σ₃_tails, σ₃_vecs=σ₃_vecs)
end
