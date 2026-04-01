@recipe(YieldSurface, C, ϕ, pT, ψ) do scene
    Attributes(
        color = :red,
        linewidth = 3,
        P_range = collect(range(-30.0, 110.0, length=1000)),
        full_dp = false
    )
end

function Makie.plot!(p::YieldSurface)
    C_obs, ϕ_obs, pT_obs = p[1], p[2], p[3]
    yield_points = lift(C_obs, ϕ_obs, pT_obs, p.full_dp) do C, ϕ, pT, full_dp
        k = sind(ϕ); c = C * cosd(ϕ); a = sqrt(1 + k^2)
        pT = full_dp ? -c / k : pT
        py = (pT + c/a) / (1 - k/a); R = py - pT; pd = py - R*k/a
        P_range = p.P_range[]; P_cap = P_range[P_range .<= pd]; P_DP = P_range[P_range .> pd]
        T_cap = sqrt.(max.(0.0, .-(P_cap .- py).^2 .+ R^2))
        T_DP = k .* P_DP .+ c
        Point2f.(vcat(P_cap, P_DP), vcat(T_cap, T_DP))
    end
    lines!(p, yield_points, color = p.color, linewidth = p.linewidth)
    return p
end
