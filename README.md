# DG MATLAB scripts
This repository contains a bunch of MATLAB code that I have been implementing and
experimenting with during the course of my MSc thesis (2019-2021) at the Delft
University of Technology.
The thesis document, titled "High-Order Discretization
of Hyperbolic Equations: Characterization of an Isogeometric Discontinuous
Galerkin Method", is publicly available [here](http://resolver.tudelft.nl/uuid:a013de4a-a869-47ff-a020-98ff5da743cd).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4767431.svg)](https://doi.org/10.5281/zenodo.4767431)

## Methods available
Three discontinous high-order discretization methods,
- Discontinuous Galerkin spectral element method (DGSEM)
- Flux reconstruction (FR)
- Modal DG with B-spline basis functions (DGIGA)

as well as optimal strong-stability-preserving explicit Runge-Kutta (SSP-RK)
time-schemes:
- 1st order, 1 stage
- 2nd order, 2 stages
- 3rd order, 3 stages
- 4th order, 5 stages
- 4th order, 10 stages

to use in combination with them.
Also, several limiters:
- Algebraic flux correction (compatible with IGA only)
- Total variation bounded slope limiter
- Biswas et al.'s moment limiter
- Burbeau et al.'s moment limiter
- Krivodonova's moment limiter
- Simple Hermite-WENO limiter

and a couple of sensors (KXRCF and Wang's AP-TVD).
On top of that, a number of support routines that extract the scheme's linear
stability and wave propagation properties (various variants of Fourier/von Neumann/modified wavenumber
analysis).

Formulation of all of these schemes is detailed in the thesis report.

## Examples
The "Examples" folder showcases the most interesting functionalities available.
The scripts in the "Experiments" directory generate all the results of the experimental part
of my master's thesis. Note that some of them might take **very** long to run.

## Installation and requirements
All you need to do is place all subdirectories of this repository in the MATLAB path
(there's a button on the "Home" tab of MATLAB's IDE that easily allows one
to do so). You will need MATLAB 2017b or newer, as well as its Parallel Computing
Toolbox.

## Author
Miquel Herrera | mail@miquelherrera.com

## License
BSD 3-clause (see `LICENSE` file)

## Cite as
[1] Miquel Herrera. mikiandh/dg-matlab-scripts: Release v1.0.0 2021. [doi:10.5281/zenodo.4767431](https://doi.org/10.5281/zenodo.4767431).

