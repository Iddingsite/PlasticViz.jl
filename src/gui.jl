function adaptive_figure_size(; frac_w = 0.72, frac_h = 0.82, min_w = 900, min_h = 700, max_w = 1600, max_h = 1200)
    try
        monitor = GLMakie.GLFW.GetPrimaryMonitor()
        mode = GLMakie.GLFW.GetVideoMode(monitor)
        sw, sh = mode.width, mode.height
        w = clamp(round(Int, sw * frac_w), min_w, max_w)
        h = clamp(round(Int, sh * frac_h), min_h, max_h)
        return (w, h)
    catch
        # Headless/unsupported backend fallback.
        return (1200, 950)
    end
end

function run_yield_plasticity(; colormap = :turbo)
    fig = Figure(size = adaptive_figure_size(), fontsize = 14)

    ui_grid = fig[2, 1:2]
    sg = SliderGrid(ui_grid[1, 1],
        (label = "Cohesion C [MPa]", range = 1.0:1.0:50.0, startvalue = 25.0),
        (label = "Friction ϕ [°]", range = 0.0:1.0:45.0, startvalue = 30.0),
        (label = "Tensile limit pT [MPa]", range = -20.0:1.0:0.0, startvalue = -5.0),
        (label = "Dilation ψ [°]", range = 0.0:1.0:45.0, startvalue = 10.0),
        width = 700, tellwidth = false
    )
    C_obs, ϕ_obs, pT_obs, ψ_obs = [s.value for s in sg.sliders]

    visible_obs = Observable(false)
    button = Button(ui_grid[2, 1][1, 1], label = @lift($visible_obs ? "Hide Construction" : "Show Construction"), width=200)
    on(button.clicks) do _
        visible_obs[] = !visible_obs[]
    end

    dp_toggle = Toggle(ui_grid[2, 1][1, 2], active = false)
    Label(ui_grid[2, 1][1, 3], "Full Drucker-Prager", fontsize = 14)

    pT_eff_obs = lift(C_obs, ϕ_obs, pT_obs, dp_toggle.active) do C, ϕ, pT, full_dp
        full_dp ? -(C * cosd(ϕ)) / sind(ϕ) : pT
    end

    orig_color_active        = sg.sliders[3].color_active[]
    orig_color_active_dimmed = sg.sliders[3].color_active_dimmed[]
    orig_color_inactive      = sg.sliders[3].color_inactive[]

    on(dp_toggle.active; update=true) do full_dp
        if full_dp
            sg.sliders[3].color_active[]        = RGBAf(0.2, 0.5, 1.0, 0.25)
            sg.sliders[3].color_active_dimmed[] = RGBAf(0.5, 0.7, 1.0, 0.25)
            sg.sliders[3].color_inactive[]      = RGBAf(0.75, 0.75, 0.75, 0.25)
            target = round(Int, pT_eff_obs[])
            sg.sliders[3].range[] = min(target, -20):1:0
            set_close_to!(sg.sliders[3], target)
        else
            sg.sliders[3].color_active[]        = orig_color_active
            sg.sliders[3].color_active_dimmed[] = orig_color_active_dimmed
            sg.sliders[3].color_inactive[]      = orig_color_inactive
            sg.sliders[3].range[] = -20.0:1.0:0.0
        end
    end

    # Update slider position when C or ϕ change while locked
    on(pT_eff_obs) do val
        dp_toggle.active[] || return
        target = round(Int, val)
        sg.sliders[3].range[] = min(target, -20):1:0
        set_close_to!(sg.sliders[3], target)
    end

    # Clamp dilation to [0, ϕ] whenever friction changes
    on(ϕ_obs; update=true) do ϕ
        if ψ_obs[] > ϕ
            set_close_to!(sg.sliders[4], floor(ϕ))
        end
    end

    # Snap dilation back if user drags beyond friction angle
    ψ_resetting = Ref(false)
    on(sg.sliders[4].value) do _
        ψ_resetting[] && return
        if ψ_obs[] > ϕ_obs[]
            ψ_resetting[] = true
            set_close_to!(sg.sliders[4], floor(ϕ_obs[]))
            ψ_resetting[] = false
        end
    end

    # Snap back if user tries to drag while locked
    pT_resetting = Ref(false)
    on(sg.sliders[3].value) do _
        dp_toggle.active[] && !pT_resetting[] || return
        pT_resetting[] = true
        set_close_to!(sg.sliders[3], round(Int, pT_eff_obs[]))
        pT_resetting[] = false
    end

    ax = Axis(fig[1, 1],
        title = "Meridional plot of the smooth yield function and flow potential",
        xlabel = "Mean Stress P [MPa]", ylabel = "Shear Stress τ [MPa]",
        aspect = DataAspect()
    )

    x_grid = range(-30.0, 110.0, length=250)
    y_grid = range(0.0, 80.0, length=250)

    # --- Q Field ---
    Q_data = lift(C_obs, ϕ_obs, pT_eff_obs, ψ_obs) do C, ϕ, pT, ψ
        k = sind(ϕ); kf = sind(ψ); c_val = C * cosd(ϕ); a = sqrt(1 + k^2); b = sqrt(1 + kf^2)
        py = (pT + c_val/a) / (1 - k/a); R = py - pT; pd = py - R*k/a; sd = k*pd + c_val
        pf = pd + kf*sd; Rf = pf - pT
        v_vec = [pd - pf, sd]; v_vec = (v_vec / (norm(v_vec) + 1e-9)) * Rf
        pdf, sdf = v_vec[1] + pf, v_vec[2]
        tol = 1e-7
        has_yield_cap = R > tol && sd > tol && abs(py - pd) > tol
        has_potential_cap = Rf > tol && sd > tol && abs(pf - pd) > tol

        [begin
            in_cap_F = has_yield_cap && (py_val * (py - pd) < (py - px) * sd)
            F = in_cap_F ? (sqrt(py_val^2 + (px - py)^2) - R) * a : (py_val - k*px - c_val)
            if F < 0.0
                NaN
            else
                in_cap_Q = has_potential_cap && (py_val * (pf - pd) < (pf - px) * sd)
                in_cap_Q ? (sqrt(py_val^2 + (px - pf)^2) - Rf) * b : (py_val - kf*(px - pdf) - sdf)
            end
        end for px in x_grid, py_val in y_grid]
    end

    cf = contourf!(ax, x_grid, y_grid, Q_data, colormap = colormap, levels = 25, nan_color = :white)

    # --- Drucker-Prager Reference ---
    dp_line = lift(C_obs, ϕ_obs) do C, ϕ
        Point2f.(x_grid, max.(0.0, sind(ϕ) .* x_grid .+ C * cosd(ϕ)))
    end
    p_dp = lines!(ax, dp_line, color = :black, linestyle = :dot, linewidth = 1.5, alpha = 0.6)

    # --- Yield Surface ---
    p_yield = yieldsurface!(ax, C_obs, ϕ_obs, pT_obs, ψ_obs, color = :red, linewidth = 4, full_dp = dp_toggle.active)

    # --- Geometric Points & Construction ---
    points_obs = lift(C_obs, ϕ_obs, pT_eff_obs, ψ_obs) do C, ϕ, pT, ψ
        k = sind(ϕ); kf = sind(ψ); c_val = C * cosd(ϕ); a = sqrt(1 + k^2)
        py = (pT + c_val/a) / (1 - k/a); R = py - pT; pd = py - R*k/a; sd = k*pd + c_val
        pf = pd + kf*sd
        return [Point2f(pT, 0), Point2f(py, 0), Point2f(pd, sd), Point2f(pf, 0)]
    end

    scatter!(ax, points_obs, color = :blue, markersize = 14, strokecolor = :white, strokewidth = 2, visible = visible_obs)
    lines!(ax, lift(p -> [p[3], p[2]], points_obs), color = :red, linestyle = :dash, linewidth = 2, visible = visible_obs)
    lines!(ax, lift(p -> [p[3], p[4]], points_obs), color = :blue, linestyle = :dash, linewidth = 2, visible = visible_obs)

    # --- Return Arrows ---
    arrow_data = lift(C_obs, ϕ_obs, pT_eff_obs, ψ_obs) do C, ϕ, pT, ψ
        p_samples = [-10.0, 10.0, 40.0, 75.0]
        k = sind(ϕ); kf = sind(ψ); c_val = C * cosd(ϕ); a = sqrt(1 + k^2)
        py = (pT + c_val/a) / (1 - k/a); R = py - pT; pd = py - R*k/a; sd = k*pd + c_val
        pf = pd + kf*sd
        tol = 1e-7
        has_yield_cap = R > tol && sd > tol && abs(py - pd) > tol
        has_potential_cap = (pf - pT) > tol && sd > tol && abs(pf - pd) > tol
        tips = Point2f[]; dirs = Vec2f[]
        for px in p_samples
            in_yield_cap = has_yield_cap && (px < pd)
            tau = in_yield_cap ? sqrt(max(0.0, R^2 - (px - py)^2)) : k*px + c_val
            push!(tips, Point2f(px, tau))
            in_potential_cap = has_potential_cap && (tau * (pf - pd) < (pf - px) * sd)
            grad = in_potential_cap ? Vec2f(px - pf, tau) : Vec2f(-kf, 1.0)
            push!(dirs, -(grad / (norm(grad) + 1e-9)) * 12.0)
        end
        return tips, dirs
    end

    arrows2d!(ax.scene, lift(x->x[1], arrow_data), lift(x->x[2], arrow_data), color = :black, align = :tip, shaftwidth = 2, tipwidth = 6)

    # --- Colorbar & Legend ---
    Colorbar(fig[1, 2], cf, label = "Plastic Potential Q [MPa]", width = 24, tellwidth = true)
    colsize!(fig.layout, 1, Auto(1))
    colsize!(fig.layout, 2, Fixed(90))
    colgap!(fig.layout, 1, 8)
    axislegend(ax,
        [p_yield, p_dp],
        ["Yield Surface", "Drucker-Prager"],
        position = :rb,
        backgroundcolor = (:white, 0.8),
        framevisible = true
    )

    hlines!(ax, 0, color = :black, alpha = 0.2); vlines!(ax, 0, color = :black, alpha = 0.2)
    xlims!(ax, -30, 110); ylims!(ax, 0, 80)
    display(fig)
end
