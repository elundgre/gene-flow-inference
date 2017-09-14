# gene-flow-inference
description a work in progress

This project is for estimating the reverse time gene flow (g) and coalescence rates (gamma)
from pairwise mean divergence and sample locations.
Currently, a grid resolution must be specified.

Scripts:

test.R
-A simple test script to demonstrate using the "run.mcmc" function 
to use the coalescence and commute time based methods to recover an example problem.
User inputs the grid size (length and width).
Gene flows rates currently are symmetric and randomly generated, but can be specified.

many.coalvcom.sym.R
-Script to run many symmetric 3x2 example problems under different situations.
Currently tests 25 different problems at 4 different noise levels
(for error in the input H).
Currently very messy

Functions:

mult_small
- a function to help parallelize the "run.mcmc" function

run.mcmc
- a function to call either the coalescence or commute time based MCMC functions three times,
once for "pre-burn-in" where the likelihood function is not as strict (to help avoid local optima),
once for normal burn-in, then once for sampling.

findG.MH
- a function that runs the Metropolis Hastings algorithm using the coalescence time method
to estimate the posterior distributions of g and gamma

findG.MH.com
- a function that runs the Metropolis Hastings algorithm using the commute time approximation
of coalescence time to estimate the posterior distributions of g and q (the within location diversity rates).

findG.nnls
- a function that uses the nnls (Non-Negative Least Squares) package to solve for g and gamma

mcmc.fns
- auxiliary functions required for the other functions, like calculating the pairwise divergence 
and log likelihood from the current values of g and gamma. Also contains some HMC functions

findG.HMC
- a Hamiltonian Monte Carlo (HMC) method based on the coalescence time. May be out of date.

