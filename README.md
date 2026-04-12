# paper-sila-ceiling

Code to reproduce: 

Gantenberg JR, La Joie R, Heston MB, and Ackley SF. Do Amyloid Trajectories Reach a Physiologic Ceiling? Evidence from Iterative Approximation and Simulation.

## Summary

This repo can be installed as an R package using `remotes` or `pak` but is designed to be loaded in place with the `pkgload` package. The package's name is _silaceiling_.

The directories are organized as follows:

- `data` summary data loaded when the package is activated, primarily simulation input parameters,
- `inst` contains the analysis scripts, including the code that produces the package datasets,
- `R` contains package assets and function definitions,
- and `output` contains analysis outputs, including simulation datasets, models fit on simulated data, and figure/table output.

Scripts in `inst` use `pkgload::load_all()` to load the package namespace. For help regarding data and functions, use a command like `?silaceiling::gen_exponential` rather than `?gen_exponential`.



