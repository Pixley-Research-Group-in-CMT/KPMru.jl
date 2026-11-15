# KPMru

Kernel Polynomial Method (KPM) calculations for tight-binding Hamiltonians: density
of states, DC conductivity, and linear and nonlinear optical response, with
optional CUDA acceleration.

[![Build Status](https://github.com/Pixley-Research-Group-in-CMT/KPMru.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/Pixley-Research-Group-in-CMT/KPMru.jl/actions/workflows/ci.yml)
[![Coverage](https://codecov.io/gh/Pixley-Research-Group-in-CMT/KPMru.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Pixley-Research-Group-in-CMT/KPMru.jl)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://pixley-research-group-in-cmt.github.io/KPMru.jl/)
[![Julia](https://img.shields.io/badge/julia-1.10%2B-blue.svg)](https://julialang.org)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22773562.svg)](https://doi.org/10.5281/zenodo.22773562)

## Installation

KPMru supports Julia 1.10 and later. In the Julia REPL, press `]` to enter package
mode, then choose one of these installation options:

**From GitHub (available now):**

```text
pkg> add https://github.com/Pixley-Research-Group-in-CMT/KPMru.jl
```

**From the General registry (registration pending):**

[Registration in General](https://github.com/JuliaRegistries/General/pull/175996)
is pending. Once it is merged and your local registry is updated, install by name:

```text
pkg> add KPMru
```

To update an existing installation, use `pkg> update KPMru`. Press Backspace at an
empty package prompt to return to Julia mode. CPU usage does not require CUDA.jl.

## Getting started with DOS

The Julia examples below run in order in the same session. They use small matrices
and expansion orders for a quick demonstration; increase the system size, expansion
order `NC`, and number of random vectors `NR` for converged calculations.
KPMru does not export its functions, so calls below use the `KPMru.` prefix.

Build a periodic nearest-neighbor chain with hopping amplitude one. The simplest
DOS interface accepts the original Hamiltonian and rescales it internally:

```julia
using KPMru
using SparseArrays
using Random

Random.seed!(1234)
NH = 16
H_orig = spdiagm(-1 => fill(-1.0, NH - 1), 1 => fill(-1.0, NH - 1))
H_orig[1, NH] = H_orig[NH, 1] = -1.0
NC = 32
NR = 4
E, rho = KPMru.dos(H_orig; NC=NC, NR=NR, N_tilde=50,
                  E_range=[-1.5, 1.5])
```

For more control, compute moments first and reuse them. `normalizeH` returns the
energy scale `a` and `H_norm = H_orig / a`, whose eigenvalues lie inside `(-1, 1)`.
The default `kpm_1d` call generates matching left and right random input vectors
for the stochastic trace estimate. Functions ending in `!` are in-place variants.

```julia
a, H_norm = KPMru.normalizeH(H_orig)
mu = KPMru.kpm_1d(H_norm, NC, NR)
rho_0 = KPMru.dos0(mu, a)
d2rho_0 = KPMru.dos0(mu, a; dE_order=2)
E_grid, rho_E = KPMru.dos(mu, a; N_tilde=50, E_range=[-1.5, 1.5])
```

`rho_0` and `d2rho_0` are the DOS and its second energy derivative at zero energy.
`N_tilde=50` divides the requested interval into 50 steps, returning 51 points
including both endpoints. `E_grid` is in the original Hamiltonian's energy units.
See [the KPM review](https://doi.org/10.1103/RevModPhys.78.275) for background.

## Response calculations

The following examples reuse the chain and parameters above. They demonstrate the
API; a small one-dimensional chain is not a model for a general response tensor.
Current operators and physical prefactors must match the model and units in an
actual calculation.

For longitudinal DC response, construct the chain current operator (in units where
the charge, lattice spacing, and reduced Planck constant are one), calculate 2D
moments, and evaluate the differential conductivity at the specified energies:

```julia
Jx = spdiagm(-1 => fill(-1.0im, NH - 1), 1 => fill(1.0im, NH - 1))
Jx[NH, 1] = 1.0im
Jx[1, NH] = -1.0im
mu_2d_xx = KPMru.kpm_2d(H_norm, Jx, Jx, NC, NR, NH)
Evals = [-0.5, 0.0, 0.5]
dsigma_xx = KPMru.d_dc_cond(mu_2d_xx, a, Evals)
```

For nonlinear response, use a smaller expansion order because the moment tensor
requires cubic storage. Frequencies and the evaluation energy passed to `d_cpge`
are in units of the rescaled Hamiltonian:

```julia
NC3 = 8
mu_3d_xxx = KPMru.kpm_3d(H_norm, Jx, Jx, Jx, NC3, NR, NH)
w1, w2, energy = 0.2, -0.15, 0.1
dchi_xxx = KPMru.d_cpge(mu_3d_xxx, NC3, w1, w2, energy)
```

These functions return differential response kernels. Obtaining a physical
conductivity requires the appropriate energy integration, occupation factors,
rescaling, and physical prefactors. See the
[documentation](https://pixley-research-group-in-cmt.github.io/KPMru.jl/),
[DC response reference](https://arxiv.org/abs/1410.8140), and
[nonlinear response references](https://arxiv.org/abs/1810.03732)
and [arXiv:2312.14244](https://arxiv.org/abs/2312.14244).

## Optional GPU support

CUDA.jl 5 is an optional dependency. GPU support is experimental and requires a
compatible NVIDIA GPU. Install CUDA separately in package mode:

```text
pkg> add CUDA@5
```

In a GPU session, load CUDA and KPMru (the following optional example is separate
from the CPU examples above):

```text
julia> using CUDA
julia> using KPMru
julia> KPMru.whichcore()
```

To query whether GPU support is active in the current session:

```julia
KPMru.whichcore()
```

Without CUDA loaded, this returns `false` and calculations run on the CPU.

## Previous package name

Earlier unregistered versions used the name `KPM`. Load this package with
`using KPMru` and update qualified calls from `KPM.function_name(...)` to
`KPMru.function_name(...)`. The moment function is named `kpm_1d` (lowercase).
