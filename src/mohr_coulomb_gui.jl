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
    secφ = sqrt(1.0 + tanφ^2)
    d    = (c + P * tanφ) / secφ          # perp. distance centre → line
    tol  = 0.5                             # half a slider step — makes :critical reachable
    state = R > d + tol ? :failed : (R > d - tol ? :critical : :safe)
    θ     = 45.0 + φ / 2.0
    # Tangent point = centre + R * unit_normal_toward_failure_line
    # Line: -tanφ·σ + τ - c = 0  →  gradient (-tanφ, 1), pointing from sub- to super-yield
    tangent_pt = Point2f(P - R * tanφ / secφ, R / secφ)
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

    ys = (-0.4f0, 0f0, 0.4f0)
    σ₃_tails = vcat([Point2f( 1f0 + σ₃_len, y) for y in ys],
                    [Point2f(-1f0 - σ₃_len, y) for y in ys])
    σ₃_vecs  = vcat([Vec2f(-σ₃_len, 0f0) for _ in ys],
                    [Vec2f( σ₃_len, 0f0) for _ in ys])

    return (plane_pts=(p1, p2), θ=θ,
            σ₁_tails=σ₁_tails, σ₁_vecs=σ₁_vecs,
            σ₃_tails=σ₃_tails, σ₃_vecs=σ₃_vecs,
            show_σ₃=σ₃ >= 0.5)
end

function run_mohr_coulomb(;
    c_default      = 15.0,
    phi_default    = 30.0,
    sigma1_default = 80.0,
    sigma3_default = 10.0
)
    fig = Figure(size = adaptive_figure_size(), fontsize = 14)

    # Row 2: controls (spans both columns)
    ui_grid = fig[2, 1:2]
    sg = SliderGrid(ui_grid[1, 1],
        (label = "Cohesion c [MPa]",      range = 0.0:1.0:50.0,  startvalue = c_default),
        (label = "Friction angle φ [°]",   range = 0.0:1.0:45.0,  startvalue = phi_default),
        (label = "Major principal stress σ₁ [MPa]",  range = 0.0:1.0:150.0, startvalue = sigma1_default),
        (label = "Minor principal stress σ₃ [MPa]",  range = 0.0:1.0:100.0, startvalue = sigma3_default),
        tellwidth = true
    )
    c_obs, φ_obs, σ₁_obs, σ₃_obs = [s.value for s in sg.sliders]

    # σ₃ ≤ σ₁ constraint (snap-back pattern from gui.jl)
    σ₃_resetting = Ref(false)
    on(sg.sliders[4].value) do _
        σ₃_resetting[] && return
        if σ₃_obs[] > σ₁_obs[]
            σ₃_resetting[] = true
            set_close_to!(sg.sliders[4], floor(σ₁_obs[]))
            σ₃_resetting[] = false
        end
    end
    on(σ₁_obs; update = true) do σ₁
        σ₃_obs[] > σ₁ || return
        σ₃_resetting[] = true
        set_close_to!(sg.sliders[4], floor(σ₁))
        σ₃_resetting[] = false
    end

    # Status label
    status_text_obs  = Observable("SAFE — stress state is inside the failure envelope")
    status_color_obs = Observable{Symbol}(:seagreen)
    Label(ui_grid[2, 1], status_text_obs,
        color = status_color_obs, fontsize = 14, font = :bold,
        tellwidth = false)
    Label(ui_grid[3, 1],
        "σ₃ is automatically limited to σ₁ — minor principal stress cannot exceed major principal stress",
        fontsize = 11, color = (:black, 0.5), font = :italic, tellwidth = false)

    # Derived observables
    mohr_state_obs   = lift(c_obs, φ_obs, σ₁_obs, σ₃_obs) do c, φ, σ₁, σ₃
        mohr_failure_state(c, φ, σ₁, σ₃)
    end
    state_sym_obs    = lift(s -> s[1], mohr_state_obs)
    tangent_obs      = lift(s -> [s[3]], mohr_state_obs)
    circle_color_obs = lift(state_sym_obs) do s
        s == :safe ? :seagreen : s == :critical ? :darkorange : :crimson
    end
    not_safe_obs = lift(s -> s != :safe, state_sym_obs)

    # Row 1: Mohr circle axis
    ax_mohr = Axis(fig[1, 1],
        title  = "Mohr Circle and Coulomb Failure Envelope",
        xlabel = "Normal stress σ [MPa]",
        ylabel = "Shear stress τ [MPa]",
        aspect = DataAspect()
    )
    hlines!(ax_mohr, 0, color = (:black, 0.25))
    vlines!(ax_mohr, 0, color = (:black, 0.25))

    # Mohr circle (parametric curve, 300 points)
    t_rng = range(0, 2π, length = 300)
    circle_pts_obs = lift(σ₁_obs, σ₃_obs) do σ₁, σ₃
        P = (σ₁ + σ₃) / 2.0
        R = (σ₁ - σ₃) / 2.0
        [Point2f(P + R * cos(t), R * sin(t)) for t in t_rng]
    end
    lines!(ax_mohr, circle_pts_obs, color = circle_color_obs, linewidth = 3)

    # Coulomb failure lines (upper and mirror below σ-axis)
    σ_ext = range(-10.0, 180.0, length = 2)  # covers max axis extent (σ₁=150 → xmax≈175)
    lines!(ax_mohr,
        lift(c_obs, φ_obs) do c, φ; [Point2f(s, c + s * tand(φ)) for s in σ_ext]; end,
        color = :black, linewidth = 2)
    lines!(ax_mohr,
        lift(c_obs, φ_obs) do c, φ; [Point2f(s, -(c + s * tand(φ))) for s in σ_ext]; end,
        color = :black, linewidth = 2, linestyle = :dash)

    # Tangent point (filled dot, visible only when not safe)
    scatter!(ax_mohr, tangent_obs,
        color = :navy, markersize = 12,
        strokecolor = :white, strokewidth = 1.5,
        visible = not_safe_obs)

    # Dashed radius from circle centre to tangent point (visible when not safe)
    lines!(ax_mohr,
        lift(σ₁_obs, σ₃_obs, mohr_state_obs) do σ₁, σ₃, s
            [Point2f((σ₁ + σ₃) / 2.0, 0), s[3]]
        end,
        color = :navy, linewidth = 1.5, linestyle = :dash,
        visible = not_safe_obs)

    # θ = 45° + φ/2 label near tangent point (visible when not safe)
    text!(ax_mohr,
        lift(s -> "θ = 45° + φ/2 = $(round(Int, s[4]))°", mohr_state_obs),
        position = lift(s -> s[3] + Vec2f(3, 3), mohr_state_obs),
        fontsize = 12, color = :navy, visible = not_safe_obs)

    # Cohesion intercept: dot at (0, c) + text label
    scatter!(ax_mohr,
        lift(c -> [Point2f(0, c)], c_obs),
        color = :darkblue, markersize = 9)
    text!(ax_mohr,
        lift(c -> "c = $(round(Int, c)) MPa", c_obs),
        position = lift(c -> Point2f(2, c + 2), c_obs),
        fontsize = 12, color = :darkblue, align = (:left, :bottom))

    # φ label on the failure line
    text!(ax_mohr,
        lift(φ -> "φ = $(round(Int, φ))°", φ_obs),
        position = lift((c, φ, σ₁) -> Point2f(σ₁ * 0.15, c + σ₁ * 0.15 * tand(φ) + 4),
                        c_obs, φ_obs, σ₁_obs),
        fontsize = 12, color = :black)

    # Failure line formula label
    text!(ax_mohr,
        lift((c, φ, σ₁) -> "τ = c + σ tan(φ)", c_obs, φ_obs, σ₁_obs),
        position = lift((c, φ, σ₁) -> Point2f(σ₁ * 0.5, c + σ₁ * 0.5 * tand(φ) + 4),
                        c_obs, φ_obs, σ₁_obs),
        fontsize = 11, color = (:black, 0.55))

    # σ₁ and σ₃ dotted vertical reference lines + axis labels
    vlines!(ax_mohr,
        lift(σ₁_obs, σ₃_obs) do σ₁, σ₃; [σ₃, σ₁]; end,
        color = (:gray, 0.4), linestyle = :dot, linewidth = 1.5)
    text!(ax_mohr, lift(σ₁ -> "σ₁=$(round(Int, σ₁))", σ₁_obs),
        position = lift(σ₁ -> Point2f(σ₁, -2), σ₁_obs),
        fontsize = 11, align = (:center, :top), color = :gray)
    text!(ax_mohr, lift(σ₃ -> "σ₃=$(round(Int, σ₃))", σ₃_obs),
        position = lift(σ₃ -> Point2f(σ₃, -2), σ₃_obs),
        fontsize = 11, align = (:center, :top), color = :gray)

    # Reactive axis limits — also account for the envelope formula label position
    function _update_mohr_limits!()
        σ₁, σ₃, c, φ = σ₁_obs[], σ₃_obs[], c_obs[], φ_obs[]
        R       = (σ₁ - σ₃) / 2.0
        y_label = c + σ₁ * 0.5 * tand(φ) + 12.0   # formula label y + margin
        xlims!(ax_mohr, min(0.0, σ₃) - 5.0, σ₁ * 1.1 + 10.0)
        ylims!(ax_mohr, -(R * 1.5 + 10.0), max(R * 1.5 + 10.0, y_label))
    end
    on(σ₁_obs) do _; _update_mohr_limits!(); end
    on(σ₃_obs) do _; _update_mohr_limits!(); end
    on(c_obs)  do _; _update_mohr_limits!(); end
    on(φ_obs)  do _; _update_mohr_limits!(); end
    _update_mohr_limits!()

    # Row 1 col 2: failure-plane block sketch
    ax_block = Axis(fig[1, 2], aspect = DataAspect())
    hidedecorations!(ax_block)
    hidespines!(ax_block)

    block_geo_obs = lift(φ_obs, σ₁_obs, σ₃_obs) do φ, σ₁, σ₃
        failure_block_geometry(φ, σ₁, σ₃)
    end

    # Material sample rectangle
    poly!(ax_block,
        [Point2f(-1,-1), Point2f(1,-1), Point2f(1,1), Point2f(-1,1)],
        color = (:steelblue, 0.12), strokecolor = :black, strokewidth = 2)

    # Failure plane (dashed red line, angle θ)
    lines!(ax_block,
        lift(geo -> [geo.plane_pts[1], geo.plane_pts[2]], block_geo_obs),
        color = :red, linewidth = 2.5, linestyle = :dash)

    # σ₁ arrows (vertical, compressive — navy)
    arrows2d!(ax_block.scene,
        lift(geo -> geo.σ₁_tails, block_geo_obs),
        lift(geo -> geo.σ₁_vecs,  block_geo_obs),
        color = :navy, shaftwidth = 2, tipwidth = 8)

    # σ₃ arrows (horizontal, compressive — orange; hidden via visible when σ₃ ≈ 0)
    σ₃_visible_obs = lift(geo -> geo.show_σ₃, block_geo_obs)
    arrows2d!(ax_block.scene,
        lift(geo -> geo.σ₃_tails, block_geo_obs),
        lift(geo -> geo.σ₃_vecs,  block_geo_obs),
        color = :darkorange, shaftwidth = 2, tipwidth = 8,
        visible = σ₃_visible_obs)

    # Stress value labels — positioned just above/right of arrow tails (reactive)
    text!(ax_block,
        lift(σ₁ -> "σ₁ = $(round(Int, σ₁)) MPa", σ₁_obs),
        position = lift(σ₁_obs) do σ₁
            l = clamp(0.15f0 + 0.4f0 * Float32(σ₁) / 150f0, 0.15f0, 0.55f0)
            Point2f(0, 1f0 + l + 0.2f0)
        end,
        align = (:center, :bottom), fontsize = 13, color = :navy)
    text!(ax_block,
        lift(σ₃ -> "σ₃ = $(round(Int, σ₃)) MPa", σ₃_obs),
        position = lift(σ₃_obs) do σ₃
            l = 0.15f0 + 0.4f0 * Float32(σ₃) / 150f0
            Point2f(1f0 + l + 0.15f0, 0)
        end,
        align = (:left, :center), fontsize = 13, color = :darkorange,
        visible = σ₃_visible_obs)

    # θ angle label near bottom-left of block
    text!(ax_block,
        lift(geo -> "θ = 45° + φ/2 = $(round(Int, geo.θ))°", block_geo_obs),
        position = Point2f(-0.9, -0.85),
        fontsize = 12, color = :red)

    # Teaching sentence below the block
    text!(ax_block,
        "The failure plane orientation depends only on φ, not on cohesion c.",
        position = Point2f(0, -1.6), align = (:center, :top),
        fontsize = 12, color = (:black, 0.65), font = :italic)

    # Reactive limits — ensure arrow tails + labels + teaching text always fit
    function _update_block_limits!()
        σ₁, σ₃ = σ₁_obs[], σ₃_obs[]
        l₁ = clamp(0.15f0 + 0.4f0 * Float32(σ₁) / 150f0, 0.15f0, 0.55f0)
        l₃ = 0.15f0 + 0.4f0 * Float32(σ₃) / 150f0
        y_top   = Float64(1f0 + l₁) + 0.8   # arrow tail + label height + margin
        x_right = Float64(1f0 + l₃) + 3.0   # arrow tail + "σ₃ = 100 MPa" text + margin
        xlims!(ax_block, -2.8, x_right)
        ylims!(ax_block, -3.0, y_top)
    end
    on(σ₁_obs) do _; _update_block_limits!(); end
    on(σ₃_obs) do _; _update_block_limits!(); end
    _update_block_limits!()

    # Wire status label to state
    on(state_sym_obs; update = true) do s
        if s == :safe
            status_text_obs[]  = "SAFE — stress state is inside the failure envelope"
            status_color_obs[] = :seagreen
        elseif s == :critical
            status_text_obs[]  = "CRITICAL — Mohr circle is tangent to the failure envelope"
            status_color_obs[] = :darkorange
        else
            status_text_obs[]  = "FAILED — Mohr circle crosses the failure envelope"
            status_color_obs[] = :crimson
        end
    end

    rowsize!(fig.layout, 1, Relative(0.82))
    rowgap!(fig.layout, 1, 8)

    display(fig)
end
