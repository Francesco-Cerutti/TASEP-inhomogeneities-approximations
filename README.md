# TASEP with Local Inhomogeneities: Approximations for the Stationary State

Comparative study of approximation methods (Mean Field, Pair Approximation, Triplet Approximation, ISA) for the stationary state of the TASEP (Totally Asymmetric Simple Exclusion Process) with local inhomogeneities (defect bonds), validated against Kinetic Monte Carlo simulations. Developed as part of my Master's thesis in Physics of Complex Systems at Politecnico di Torino.

## Model

The TASEP describes particles hopping unidirectionally along a 1D lattice of `L` sites, entering at rate `α`, exiting at rate `β`, and hopping between internal sites at rate `q` (1 for a homogeneous chain). A **defect** is a stretch of one or more bonds with a reduced hopping rate `q < 1`, which locally slows the particle current and can induce phase transitions in the stationary density profile.

The repo implements and compares four ways of computing the stationary current/density:
- **Mean Field (MF)** — factorized approximation, fastest but least accurate near defects.
- **Pair Approximation (PA)** — accounts for nearest-neighbor correlations.
- **Triplet Approximation (TRI)** — accounts for next-nearest-neighbor correlations, most accurate of the analytical approximations.
- **ISA (Interacting Subsystem Approximation)** — closed-form approximation specific to the defect region.
- **Kinetic Monte Carlo (Gillespie algorithm)** — exact stochastic simulation, used as ground truth.

## Contents

| File | Description |
|---|---|
| `oggetti.jl` | Core library: `Tasep` struct, lattice construction with defects (`lattice`, `centered_positions`), and the Euler-integration solvers for MF, PA, and Triplet approximations, plus the ISA closed-form solution and correlation utilities. |
| `gillespie_temp.jl` | Event-driven (Gillespie) kinetic Monte Carlo simulator for the TASEP, used both as a standalone method and as the numerical benchmark for the approximations. |
| `Corrente_in_funzione_di_q.ipynb` | Current `J` as a function of the defect strength `q`, for defects of increasing length (1 to 7 bonds), comparing all methods against Kinetic Monte Carlo. |
| `Corrente_in_funzione_della_densità.ipynb` | Reproduces the fundamental diagram `J(ρ) = ρ(1-ρ)` for the pure and defected TASEP, comparing all methods. |
| `density_profiles_compared.ipynb` | Stationary density profiles `ρ(site)` across the different phases (low-density, high-density, coexistence, maximal-current), for both the pure TASEP and the TASEP with a central defect. |
| `relaxation_time_gillespie.ipynb` | Determines the thermalization (relaxation) time needed for the Gillespie simulation to reach the stationary state, across the different phases, by tracking the convergence of the current, of P10, and of the density profile. |

## Requirements

- [Julia](https://julialang.org/) (tested on 1.x)
- Packages: `OffsetArrays`, `Plots`, `LinearAlgebra`, `Statistics`, `Revise`, `ProgressMeter`, `JLD2`, `SpecialFunctions`, `Random`, `BenchmarkTools`, `Profile`, `LaTeXStrings`, `Roots`

Install the dependencies from the Julia REPL:
```julia
import Pkg
Pkg.add(["OffsetArrays", "Plots", "LinearAlgebra", "Statistics", "Revise", "ProgressMeter", "JLD2", "SpecialFunctions", "Random", "BenchmarkTools", "Profile", "LaTeXStrings", "Roots"])
```

## Usage

1. Keep `oggetti.jl` and `gillespie_temp.jl` in the same directory as the notebooks — they're loaded via `include("oggetti.jl")` and `include("gillespie_temp.jl")`.
2. Open a notebook in Jupyter (IJulia kernel) or VS Code with the Julia extension, and run the cells in order.
3. Each analysis notebook loads pre-computed results from `.jld2` data files (saved with `@save`/loaded with `@load`) rather than always re-running the simulations, since some Monte Carlo runs are expensive. To regenerate a result from scratch, uncomment the relevant simulation call and comment out the corresponding `@load` line.

## Notes

- Cell outputs are cleared before committing to keep the repository lightweight — re-run a notebook to regenerate its plots.
- Plots use a static backend (`gr()`); avoid switching to an interactive backend (e.g. `plotlyjs()`) before saving outputs, as this can produce very large notebook files.
- Saved `.jld2` result files and exported plot images are expected under a local `data/` and `images/` folder structure (see the `path_data` / `path_image` variables at the top of each notebook) — these folders are not tracked in the repo and should be created locally, or the paths adjusted to your own setup.

## Author

Francesco Cerutti — Master's student in Physics of Complex Systems, Politecnico di Torino / Università di Torino.
