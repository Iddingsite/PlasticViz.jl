# PlasticViz.jl

**PlasticViz.jl** is an interactive Julia tool for visualising yield surfaces in the meridional stress space (mean pressure P – shear stress τ). It implements the smooth linearised Drucker-Prager shear envelope combined with a circular tensile cap, following the formulation of Popov et al. (2025).

In plasticity modelling, the P–τ meridional plane is the standard diagnostic space: it shows how a material transitions from sub-yield behaviour (elastic or viscous) to plastic flow as confining pressure and deviatoric stress change. PlasticViz lets you explore this graphically in real time:

- The **yield surface** separates the sub-yield domain (inside) from the plastic domain (outside). Its shape is controlled by cohesion, friction angle, and tensile limit. Stress states inside the yield surface correspond to **elastic or viscous** behaviour, no irreversible plastic strain is produced.
- The **plastic potential field** (background colour) encodes the scalar distance to the yield surface: positive values indicate yielded material, negative values sub-yield (elastic or viscous) material.
- The **return-direction arrows** illustrate where a stress state outside the yield surface would be mapped back during a closest-point return-mapping algorithm. Their inclination relative to the yield surface normal reflects the degree of non-associativity, governed by the dilation angle.

The package is intended for teaching and research in geomechanics and computational plasticity. It is not a standalone constitutive solver.

## Features


![Example usage](media/example_usage.gif)


- Interactive sliders for cohesion (C), friction angle (ϕ), tensile limit (pT), and dilation angle (ψ).
- Optional Full Drucker-Prager mode (tensile limit locked to cone apex).
- Dropdown menu to switch colormaps interactively, with an invert toggle.
- Yield surface, reference DP line, and plastic potential field visualisation.
- Construction geometry overlay and return-direction arrows (togglable).
- **CSV import**: drag-and-drop or load a `.csv` file of stress points (P [MPa]; Shear [MPa]) to overlay them on the diagram. Points are automatically classified as:
  - *Elastic / Viscous* — inside the yield surface (no plastic yielding)
  - *Mode 1* — plastic on the tensile cap
  - *Mode 2* — plastic on the Drucker-Prager shear envelope
- Classification and plastic-point counter update live when sliders are adjusted.

## Installation

In Julia:

```julia
using Pkg
Pkg.add(url="https://github.com/Iddingsite/PlasticViz.jl")
```

## Usage

```julia
using PlasticViz

run_yield_plasticity()

# Choose a different startup colormap
run_yield_plasticity(colormap = :viridis)
```

## Parameters

| Parameter | Symbol | Unit | Description |
|-----------|--------|------|-------------|
| Cohesion | C | MPa | Shear strength of the material at zero confining pressure. It sets the intercept of the Drucker-Prager envelope on the τ axis. Higher cohesion shifts the yield surface upward, allowing the material to sustain more shear stress before yielding. |
| Friction angle | ϕ | ° | Controls the slope of the linear Drucker-Prager shear envelope. A higher friction angle means yield strength increases more steeply with confining pressure, typical of granular or frictional geomaterials. The dilation angle cannot exceed this value. |
| Tensile limit | pT | MPa | Sets the location of the tensile cut-off: the maximum tensile mean stress the material can sustain before the circular cap closes the yield surface. When Full Drucker-Prager is active, this is computed automatically as the apex of the cone (`-C cos ϕ / sin ϕ`) and the slider is locked. |
| Dilation angle | ψ | ° | Controls the direction of plastic flow (the flow potential). When ψ = ϕ the flow is associative (plastic strain normal to yield surface). When ψ < ϕ the flow is non-associative, producing less volumetric expansion during shearing, which is more realistic for most geomaterials. |

## Reference

Popov, A. A., Berlie, N., & Kaus, B. J. (2025). A dilatant visco-elasto-viscoplasticity model with globally continuous tensile cap: stable two-field mixed formulation. Geoscientific Model Development, 18(19), 7035-7058.