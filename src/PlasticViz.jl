module PlasticViz

using GLMakie
using LinearAlgebra

export run_yield_plasticity
export run_mohr_coulomb

include("recipe.jl")
include("gui.jl")
include("mohr_coulomb_gui.jl")

function __init__()
    # GLMakie opens a hidden test window during OpenGL backend initialisation.
    # Close it immediately so no blank window appears on `using PlasticViz`.
    GLMakie.closeall()
end

end # module PlasticViz
