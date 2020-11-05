# Summary
The MATLAB program provided here (bacterialtransportsimulation.m) numerically solves the coupled bacterial transport equations to model the distribution of bacterial cells exposed to a gradient of chemoattractant. Additionally the program implements a fitting procedure which can fit the model to experimental data and provide estimates for input parameters which are not obtained experimentally, specifically the random motility coefficient, the receptor/ligand dissociation constant, and the chemotactic sensitivity coefficient.

This program was produced as part of a Master of Physics project which investigated antibiotic challenge of spectinomycin and kanamycin for *Bacillus subtilis* swimming velocity and chemotaxis.

# Details
The *pdepe* solver is used to solve the bacterial transport equations, which uses a variable-step variable-order solver to perform time integration (1). A homogenous intial distribution of bacteria is used as the initial condition and no-flux Neumann boundary conditions are applied to both boundaries.
The fitting procedure minimises the discrepancy between predicted and measured bacteria distribution profiles, calculated as the sum of the squared error, using *fminsearch* which uses a Nelder-Mead simplex algorithm (2). Finally, the sum of squared error was found to be reduced significantly (roughly by a factor of 3) by disregarding the "tail" in the bacterial distribution profile from population exposed to a gradient of chemoattractant. Hence, the program also includes the option to disregard this "tail" when fitting.

# Referenced Literature
1. Shampine, L. F. & Reichelt, M. W. The MATLAB ODE suite. *SIAM J. Sci. Comput.* **18**, 1–22 (1997).
2. Lagarias, J. C., Reeds, J. A., Wright, M. H. & Wright, P. E. Convergence properties of Nelder-Mead simplex method in low dimensions. *SIAM J. Optim.* **9**, 112–147 (1998).
