# Helper: classify a Mohr-Coulomb stress state and return diagnostic geometry.
# Returns: (state::Symbol, d::Float64, tangent_pt::Point2f, θ::Float64)
#   state       — :safe | :critical | :tensile
#   d           — perpendicular distance from circle centre to failure line
#   tangent_pt  — closest point on circle to the failure line
#   θ           — failure plane angle in degrees (45 + φ/2)
function mohr_failure_state(c::Real, φ::Real, σ₁::Real, σ₃::Real; T₀::Real=Inf)
    P    = (σ₁ + σ₃) / 2.0
    R    = (σ₁ - σ₃) / 2.0
    tanφ = tand(φ)
    secφ = sqrt(1.0 + tanφ^2)
    d    = (c + P * tanφ) / secφ          # perp. distance centre → line
    tol  = 0.1                             # half a slider step — makes boundary states reachable
    state = if σ₃ ≤ -T₀ + tol
        :tensile
    elseif R > d - tol
        :critical
    else
        :safe
    end
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
    σ₃_len = max(0.05f0, 0.15f0 + 0.4f0 * Float32(σ₃) / 150f0)

    xs = (-0.5f0, 0f0, 0.5f0)
    if σ₁ >= 0
        # Compressive: arrows point into top/bottom faces
        σ₁_tails = vcat([Point2f(x,  1f0 + σ₁_len) for x in xs],
                        [Point2f(x, -1f0 - σ₁_len) for x in xs])
        σ₁_vecs  = vcat([Vec2f(0f0, -σ₁_len) for _ in xs],
                        [Vec2f(0f0,  σ₁_len) for _ in xs])
    else
        # Tensile: arrows point away from top/bottom faces
        σ₁_tails = vcat([Point2f(x,  1f0) for x in xs],
                        [Point2f(x, -1f0) for x in xs])
        σ₁_vecs  = vcat([Vec2f(0f0,  σ₁_len) for _ in xs],
                        [Vec2f(0f0, -σ₁_len) for _ in xs])
    end

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
    sigma1_default = 60.0,
    sigma3_default = 10.0
)
    _fsz = adaptive_figure_size()
    fig = Figure(size = (_fsz[2], _fsz[2]), fontsize = 12)

    # Row 2: controls (spans both columns)
    ui_grid = fig[2, 1:2]
    sg = SliderGrid(ui_grid[1, 1],
        (label = "Cohesion c [MPa]",                  range = 0.0:0.1:50.0,  startvalue = c_default),
        (label = "Friction angle φ [°]",               range = 0.0:0.1:45.0,  startvalue = phi_default),
        (label = "Major principal stress σ₁ [MPa]",   range = -30.0:0.1:150.0, startvalue = sigma1_default),
        (label = "Minor principal stress σ₃ [MPa]",   range = -30.0:0.1:100.0, startvalue = sigma3_default),
        (label = "Tensile strength T₀ [MPa]",         range = 0.0:0.1:30.0,  startvalue = 0.0),
        (label = "Pore pressure Pf [MPa]",             range = 0.0:0.1:100.0, startvalue = 0.0),
        tellwidth = true
    )
    c_obs, φ_obs, σ₁_obs, σ₃_obs, T₀_obs, u_obs = [s.value for s in sg.sliders]

    # Effective stresses (failure criterion always in effective stress space)
    σ₁eff_obs = lift((σ₁, u) -> σ₁ - u, σ₁_obs, u_obs)
    σ₃eff_obs = lift((σ₃, u) -> σ₃ - u, σ₃_obs, u_obs)

    # σ₃ two-sided constraint
    #   Ceiling — σ₃ ≤ σ₁ (total stresses)
    #   Floor   — σ₃eff ≥ −T₀: leftmost point of circle cannot cross the tensile cap
    σ₃_resetting = Ref(false)
    function _snap_σ₃!()
        σ₃_resetting[] && return
        u, T₀ = u_obs[], T₀_obs[]
        if σ₃_obs[] > σ₁_obs[]
            σ₃_resetting[] = true
            set_close_to!(sg.sliders[4], floor(σ₁_obs[] * 10 + 1e-9) / 10)
            σ₃_resetting[] = false
            return
        end
        σ₃_floor = u - T₀
        if σ₃_obs[] < σ₃_floor
            σ₃_resetting[] = true
            set_close_to!(sg.sliders[4], ceil(σ₃_floor * 10 - 1e-9) / 10)
            σ₃_resetting[] = false
        end
    end
    on(σ₃_obs) do _; _snap_σ₃!(); end
    on(σ₁_obs) do _; _snap_σ₃!(); end
    on(u_obs)  do _; _snap_σ₃!(); end
    on(T₀_obs) do _; _snap_σ₃!(); end
    _snap_σ₃!()

    # σ₁ cap: two-sided constraint
    #   Floor — σ₁eff ≥ −T₀: prevents the whole circle from sitting left of the tensile cutoff
    #   Ceiling — σ₁eff ≤ shear failure threshold (compressive regime only)
    σ₁_resetting = Ref(false)
    function _snap_σ₁!()
        σ₁_resetting[] && return
        c, φ, σ₃, u, T₀ = c_obs[], φ_obs[], σ₃_obs[], u_obs[], T₀_obs[]

        # Floor: σ₁eff = σ₁ − u ≥ −T₀  →  σ₁ ≥ u − T₀
        σ₁_floor = u - T₀
        if σ₁_obs[] < σ₁_floor
            σ₁_resetting[] = true
            set_close_to!(sg.sliders[3], ceil(σ₁_floor * 10 - 1e-9) / 10)
            σ₁_resetting[] = false
            return
        end

        # Ceiling: shear failure threshold
        s = sind(φ)
        abs(s - 1.0) < 1e-6 && return
        σ₃eff     = σ₃ - u
        σ₁eff_max = (2*c*cosd(φ) + σ₃eff*(1+s)) / (1-s)
        σ₁_max    = max(σ₃, σ₁eff_max + u, σ₁_floor)   # never push σ₁ below tensile floor
        σ₁_obs[] ≤ σ₁_max && return
        σ₁_resetting[] = true
        # Clamp to σ₃ after rounding to avoid FP edge case where floor drops below σ₃
        snap_val = max(floor(σ₁_max * 10 + 1e-9) / 10, σ₃_obs[])
        set_close_to!(sg.sliders[3], snap_val)
        σ₁_resetting[] = false
    end
    on(σ₁_obs) do _; _snap_σ₁!(); end
    on(c_obs)  do _; _snap_σ₁!(); end
    on(φ_obs)  do _; _snap_σ₁!(); end
    on(σ₃_obs) do _; _snap_σ₁!(); end
    on(u_obs)  do _; _snap_σ₁!(); end
    on(T₀_obs) do _; _snap_σ₁!(); end

    # Status label
    status_text_obs  = Observable("SAFE — stress state is inside the failure envelope")
    status_color_obs = Observable{Symbol}(:seagreen)
    Label(ui_grid[2, 1], status_text_obs,
        color = status_color_obs, fontsize = 14, font = :bold,
        tellwidth = false)
    Label(ui_grid[3, 1],
        "σ₃ ≤ σ₁ always.  σ₁ is capped at the shear failure threshold.  Failure criterion uses effective stresses σ − u.",
        fontsize = 11, color = (:black, 0.5), font = :italic, tellwidth = false)

    # Derived observables
    mohr_state_obs = lift(c_obs, φ_obs, σ₁eff_obs, σ₃eff_obs, T₀_obs) do c, φ, σ₁eff, σ₃eff, T₀
        mohr_failure_state(c, φ, σ₁eff, σ₃eff; T₀=T₀)
    end
    state_sym_obs    = lift(s -> s[1], mohr_state_obs)
    tangent_obs      = lift(s -> [s[3]], mohr_state_obs)
    circle_color_obs = lift(state_sym_obs) do s
        s == :safe ? :seagreen : s == :critical ? :red : :purple
    end
    # Tangent point decoration only meaningful for shear failure, not tensile
    shear_not_safe_obs = lift(s -> s == :critical, state_sym_obs)

    # Row 1 col 1: Mohr circle axis
    ax_mohr = Axis(fig[1, 1],
        title  = "Mohr Circle and Coulomb Failure Envelope",
        xlabel = "Normal stress σ [MPa]",
        ylabel = "Shear stress τ [MPa]",
        aspect = DataAspect(),
        backgroundcolor = RGBf(0.94, 0.96, 0.99)
    )

    text!(ax_mohr,
        0.02, 0.97,
        text = L"\tau_{\max} = \frac{\sigma_1 - \sigma_3}{2}",
        space = :relative,
        align = (:left, :top),
        fontsize = 12
    )
    hlines!(ax_mohr, 0, color = (:black, 0.25))
    vlines!(ax_mohr, 0, color = (:black, 0.25))

    # Mohr circles
    t_rng = range(0, 2π, length = 300)

    # Total stress circle — faint dotted, shown only when u > 0
    circle_pts_tot_obs = lift(σ₁_obs, σ₃_obs) do σ₁, σ₃
        P = (σ₁ + σ₃) / 2.0
        R = (σ₁ - σ₃) / 2.0
        [Point2f(P + R * cos(t), R * sin(t)) for t in t_rng]
    end
    show_total_obs = lift(u -> u > 0.5, u_obs)
    lines!(ax_mohr, circle_pts_tot_obs,
        color = (:gray, 0.35), linewidth = 1.5, linestyle = :dot,
        visible = show_total_obs)

    # Effective stress circle — main, solid, coloured by failure state
    circle_pts_eff_obs = lift(σ₁eff_obs, σ₃eff_obs) do σ₁eff, σ₃eff
        P = (σ₁eff + σ₃eff) / 2.0
        R = (σ₁eff - σ₃eff) / 2.0
        [Point2f(P + R * cos(t), R * sin(t)) for t in t_rng]
    end
    lines!(ax_mohr, circle_pts_eff_obs, color = circle_color_obs, linewidth = 3)

    # Coulomb failure envelope: smooth circular arc at tensile cutoff + straight Coulomb line.
    # Circle centred at (σ_py, 0) with radius R_cap passes through (−T₀, 0) and is tangent
    # to τ = c + σ·tan(φ) at the pivot point (σ_py − R_cap·sinφ, R_cap·cosφ).
    function _envelope_upper(c, φ, T₀)
        sinφ = sind(φ); cosφ = cosd(φ); tanφ = tand(φ)
        R_cap = (c*cosφ - T₀*sinφ) / (1.0 - sinφ)
        if R_cap > 1e-6
            σ_py  = R_cap - T₀
            θ_end = atan(cosφ, -sinφ)          # angle at pivot relative to centre
            θs    = range(Float64(π), θ_end, length = 80)
            arc   = [Point2f(σ_py + R_cap*cos(θ), R_cap*sin(θ)) for θ in θs]
            return vcat(arc, [Point2f(200.0, c + 200.0*tanφ)])
        else
            return [Point2f(-150.0, c - 150.0*tanφ), Point2f(200.0, c + 200.0*tanφ)]
        end
    end
    envelope_upper_obs = lift(_envelope_upper, c_obs, φ_obs, T₀_obs)
    envelope_lower_obs = lift(pts -> [Point2f(p[1], -p[2]) for p in pts], envelope_upper_obs)
    lines!(ax_mohr, envelope_upper_obs, color = :black, linewidth = 2)
    lines!(ax_mohr, envelope_lower_obs, color = :black, linewidth = 2, linestyle = :dash)

    # Tracks the current bottom y limit of the Mohr axis (updated by _update_mohr_limits!)
    y_bottom_obs = Observable(-50.0)

    # Tensile cutoff — vertical dashed line at σ = −T₀ (always shown)
    vlines!(ax_mohr, lift(T₀ -> [-T₀], T₀_obs),
        color = :purple, linewidth = 2, linestyle = :dash)
    text!(ax_mohr,
        lift(T₀ -> "T₀ = $(round(Int, T₀)) MPa", T₀_obs),
        position = lift((T₀, yb) -> Point2f(-T₀ + 1.5, yb + 2.0), T₀_obs, y_bottom_obs),
        fontsize = 11, color = :purple, align = (:left, :bottom))

    # Tangent point (shear failure only)
    scatter!(ax_mohr, tangent_obs,
        color = :red, markersize = 12,
        strokecolor = :white, strokewidth = 1.5,
        visible = shear_not_safe_obs)

    # Dashed radius from circle centre to tangent point (shear failure only)
    lines!(ax_mohr,
        lift(σ₁eff_obs, σ₃eff_obs, mohr_state_obs) do σ₁eff, σ₃eff, s
            [Point2f((σ₁eff + σ₃eff) / 2.0, 0), s[3]]
        end,
        color = :red, linewidth = 1.5, linestyle = :dash,
        visible = shear_not_safe_obs)

    # θ label near tangent point (shear failure only)
    text!(ax_mohr,
        lift(s -> "θ = 45° + φ/2 = $(round(Int, s[4]))°", mohr_state_obs),
        position = lift(s -> s[3] + Vec2f(3, 3), mohr_state_obs),
        fontsize = 12, color = :red, visible = shear_not_safe_obs)

    # Cohesion intercept: dot at (0, c) + text label
    scatter!(ax_mohr,
        lift(c -> [Point2f(0, c)], c_obs),
        color = :darkblue, markersize = 9)
    text!(ax_mohr,
        lift(c -> "c = $(round(Int, c)) MPa\n ", c_obs),
        position = lift(c -> Point2f(2.0, c + 3), c_obs),
        fontsize = 12, color = :darkblue, align = (:left, :bottom))

    # Failure line formula label
    text!(ax_mohr,
        lift((c, φ, σ₁) -> "τ = c + σ tan(φ)\n ", c_obs, φ_obs, σ₁_obs),
        position = lift((c, φ, σ₁) -> Point2f(σ₁ * 0.5, c + σ₁ * 0.5 * tand(φ) + 4),
                        c_obs, φ_obs, σ₁_obs),
        fontsize = 11, color = (:black, 0.55))

    # σ₁ and σ₃ effective stress reference lines + labels
    # (show σ' notation when pore pressure is active)
    vlines!(ax_mohr,
        lift(σ₁eff_obs, σ₃eff_obs) do σ₁eff, σ₃eff; [σ₃eff, σ₁eff]; end,
        color = (:gray, 0.4), linestyle = :dot, linewidth = 1.5)
    text!(ax_mohr,
        lift((σ₁eff, u) -> u > 0.5 ? "σ₁' = $(round(Int, σ₁eff)) MPa" : "σ₁ = $(round(Int, σ₁eff)) MPa",
             σ₁eff_obs, u_obs),
        position = lift(σ₁eff -> Point2f(σ₁eff + 2.0, 2.0), σ₁eff_obs),
        fontsize = 11, align = (:left, :bottom), color = :gray)
    text!(ax_mohr,
        lift((σ₃eff, u) -> u > 0.5 ? "σ₃' = $(round(Int, σ₃eff)) MPa" : "σ₃ = $(round(Int, σ₃eff)) MPa",
             σ₃eff_obs, u_obs),
        position = lift(σ₃eff -> Point2f(σ₃eff + 2.0, -2.0), σ₃eff_obs),
        fontsize = 11, align = (:left, :top), color = :gray)

    # Reactive axis limits — equal x/y range so DataAspect gives stable layout
    function _update_mohr_limits!()
        σ₁, σ₃, c, φ, u, T₀ = σ₁_obs[], σ₃_obs[], c_obs[], φ_obs[], u_obs[], T₀_obs[]
        R     = (σ₁ - σ₃) / 2.0
        σ₃eff = σ₃ - u
        x_lo  = min(0.0, σ₃eff, -T₀) - 8.0
        x_hi  = max(σ₁, σ₁ - u + R) * 1.1 + 10.0
        y_hi  = max(R * 1.5 + 10.0, c + σ₁ * 0.5 * tand(φ) + 12.0)
        x_range = x_hi - x_lo
        y_range = 2.0 * y_hi
        if x_range < y_range
            extra = (y_range - x_range) / 2.0
            x_lo -= extra;  x_hi += extra
        else
            y_hi = x_range / 2.0
        end
        xlims!(ax_mohr, x_lo, x_hi)
        ylims!(ax_mohr, -y_hi, y_hi)
        y_bottom_obs[] = -y_hi
    end
    on(σ₁_obs) do _; _update_mohr_limits!(); end
    on(σ₃_obs) do _; _update_mohr_limits!(); end
    on(c_obs)  do _; _update_mohr_limits!(); end
    on(φ_obs)  do _; _update_mohr_limits!(); end
    on(u_obs)  do _; _update_mohr_limits!(); end
    on(T₀_obs) do _; _update_mohr_limits!(); end
    _update_mohr_limits!()

    # Row 1 col 2: failure-plane block sketch (uses total stresses for applied loads)
    ax_block = Axis(fig[1, 2])
    hidedecorations!(ax_block)
    hidespines!(ax_block)

    block_geo_obs = lift(φ_obs, σ₁_obs, σ₃_obs) do φ, σ₁, σ₃
        failure_block_geometry(φ, σ₁, σ₃)
    end

    # Material sample rectangle
    poly!(ax_block,
        [Point2f(-1,-1), Point2f(1,-1), Point2f(1,1), Point2f(-1,1)],
        color = (:steelblue, 0.12), strokecolor = :black, strokewidth = 2)

    is_tensile_obs   = lift(s -> s == :tensile, state_sym_obs)
    not_tensile_obs  = lift(s -> s != :tensile, state_sym_obs)

    # Shear failure plane (diagonal dashed, faint when safe) — hidden during tensile
    plane_color_obs = lift(state_sym_obs) do s
        s == :safe ? (:red, 0.2) : (:red, 1.0)
    end
    lines!(ax_block,
        lift(geo -> [geo.plane_pts[1], geo.plane_pts[2]], block_geo_obs),
        color = plane_color_obs, linewidth = 2.5, linestyle = :dash,
        visible = not_tensile_obs)

    # Tensile crack — vertical solid crimson line (Mode I, perpendicular to σ₃)
    lines!(ax_block, [Point2f(0, -1), Point2f(0, 1)],
        color = :purple, linewidth = 3.0, linestyle = :solid,
        visible = is_tensile_obs)

    # σ₁ arrows (vertical, compressive — navy; always shown)
    arrows2d!(ax_block.scene,
        lift(geo -> geo.σ₁_tails, block_geo_obs),
        lift(geo -> geo.σ₁_vecs,  block_geo_obs),
        color = :navy, shaftwidth = 2, tipwidth = 8)

    # σ₃ compressive arrows (inward, orange) — hidden during tensile failure or σ₃ ≈ 0
    σ₃_comp_visible_obs = lift((geo, s) -> geo.show_σ₃ && s != :tensile,
                               block_geo_obs, state_sym_obs)
    arrows2d!(ax_block.scene,
        lift(geo -> geo.σ₃_tails, block_geo_obs),
        lift(geo -> geo.σ₃_vecs,  block_geo_obs),
        color = :darkorange, shaftwidth = 2, tipwidth = 8,
        visible = σ₃_comp_visible_obs)

    # σ₃ tensile arrows (outward from block faces, purple) — only during tensile failure
    # Fixed geometry: tails on block faces, vectors pointing away
    _ys_t = (-0.4f0, 0f0, 0.4f0)
    σ₃_tens_tails = vcat([Point2f( 1f0, y) for y in _ys_t],
                         [Point2f(-1f0, y) for y in _ys_t])
    σ₃_tens_vecs  = vcat([Vec2f( 0.35f0, 0f0) for _ in _ys_t],
                         [Vec2f(-0.35f0, 0f0) for _ in _ys_t])
    arrows2d!(ax_block.scene, σ₃_tens_tails, σ₃_tens_vecs,
        color = :purple, shaftwidth = 2, tipwidth = 8,
        visible = is_tensile_obs)

    # Stress value labels
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
        visible = σ₃_comp_visible_obs)

    # θ label (shear failure only) — top-left corner of block
    text!(ax_block,
        lift(geo -> "θ = 45° + φ/2 = $(round(Int, geo.θ))°", block_geo_obs),
        position = Point2f(-0.95, 0.92),
        fontsize = 12, color = :red, align = (:left, :top),
        visible = not_tensile_obs)

    # Tensile failure label (replaces θ label)
    text!(ax_block, "Mode I tensile crack\n(opens ⊥ to σ₃)",
        position = Point2f(-0.95, 0.92),
        fontsize = 12, color = :purple, align = (:left, :top),
        visible = is_tensile_obs)

    # Teaching sentence below the block
    text!(ax_block,
        "The failure plane orientation depends\n only on φ, not on cohesion c.",
        position = Point2f(0, -1.6), align = (:center, :top),
        fontsize = 12, color = (:black, 0.65), font = :italic,
        visible = not_tensile_obs)
    text!(ax_block,
        "Tensile failure: effective σ₃ falls\n below the tensile strength −T₀.",
        position = Point2f(0, -1.6), align = (:center, :top),
        fontsize = 12, color = (:purple, 0.8), font = :italic,
        visible = is_tensile_obs)

    # Reactive limits — tighter margins so the block fills more of the panel
    function _update_block_limits!()
        σ₁, σ₃ = σ₁_obs[], σ₃_obs[]
        l₁ = clamp(0.15f0 + 0.4f0 * Float32(σ₁) / 150f0, 0.15f0, 0.55f0)
        l₃ = 0.15f0 + 0.4f0 * Float32(σ₃) / 150f0
        y_top   = Float64(1f0 + l₁) + 0.6
        x_right = Float64(1f0 + l₃) + 1.8
        xlims!(ax_block, -1.8, x_right)
        ylims!(ax_block, -2.5, y_top)
    end
    on(σ₁_obs) do _; _update_block_limits!(); end
    on(σ₃_obs) do _; _update_block_limits!(); end
    _update_block_limits!()

    # Wire status label to state
    on(state_sym_obs; update = true) do s
        if s == :safe
            status_text_obs[]  = "NO FAILURE — stress state is inside the failure envelope"
            status_color_obs[] = :seagreen
        elseif s == :critical
            status_text_obs[]  = "FAILURE — Mohr circle is tangent to the failure envelope"
            status_color_obs[] = :red
        elseif s == :tensile
            status_text_obs[]  = "TENSILE FAILURE — effective minor stress exceeds tensile strength T₀"
            status_color_obs[] = :purple
        end
    end

    rowsize!(fig.layout, 1, Relative(0.68))
    colsize!(fig.layout, 1, Relative(0.45))
    rowgap!(fig.layout, 1, 8)

    display(fig)
end
