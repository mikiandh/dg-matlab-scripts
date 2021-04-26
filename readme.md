# DG matlab scripts
This repository contains a bunch of Matlab code that I have used during my MSc graduation project (2019-2021), titled: "High-order discretization of hyperbolic PDEs: characterization of an isogeometric discontinuous Galerkin method".

This includes three discontinous high-order discretization methods,
1. Discontinuous Galerkin spectral element method (DGSEM)
2. Flux reconstruction (FR)
3. Modal DG with B-spline basis functions (DGIGA)

as well as explicit Runge-Kutta time-schemes to use in combination with them.
Also, a number of support routines I used to extract their linear stability and wave propagation properties, several limiters
1. Algebraic flux correction (compatible with IGA only)
2. Total variation bounded slope limiter
3. Biswas et al.'s moment limiter
4. Burbeau et al.'s moment limiter
5. Krivodonova's moment limiter
6. Simple Hermite-WENO limiter

and a couple of sensors (KXRCF and Wang's AP-TVD).

## Examples
The "Examples" folder showcases the most interesting functionalities implemented.
The scripts in the "Experiments" directory generate all the results of the experimental part of my master's thesis. Note that some of them might take **very** long to run.

## Installation
You will need MATLAB 2017b or newer.
All you need to do is place all subdirectories of this repository in the MATLAB path (there's a button on the "Home" tab that easily allows one to do so).

## Author
Miquel Herrera | mail@miquelherrera.com

## License
BSD 3-clause (see `LICENSE` file)