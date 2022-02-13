# BinaryStars
[![License: New BSD](https://img.shields.io/badge/License-NewBSD-yellow.svg)](https://opensource.org/licenses/BSD-3-Clause)

An efficient Python implementation for Bayesian inference in binary stars based on the probabilistic programming language [Stan](https://mc-stan.org/).

![alt text](https://github.com/mvidela31/BinaryStars/blob/main/examples/results/HIP109951_observations_plot.png)
![alt text](https://github.com/mvidela31/BinaryStars/blob/main/examples/results/HIP109951_params_plot.png)

## Contents
This repository includes the implementation of the statistical models for different binary stellar systems configurations:
- Visual Binary model (`visual.stan`).
- SB2 model (`sb2.stan`).
- Visual-SB2 model (`visual_sb2.stan`).
- Visual-SB1 model (`visual_sb1.stan`).

## Requirements
* [CmdStanPy](https://mc-stan.org/cmdstanpy/).
* [Python 3](https://www.python.org/).

## Usage
See the `examples/Visual_SB2.ipynb` notebook for the instructions usage. To incorporate prior distributions see the `examples/Visual_SB1_priors.ipynb` notebook example.

## References
[1]  Videla, M., Mendez, R. A., Claveria, R. M., Silva, J. F., & Orchard, M. E. (in prep). **Bayesian inference in single-line spectroscopic binaries with a visual orbit**.

[[2](https://www.osti.gov/biblio/1430202)] Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M., ... & Riddell, A. (2017). **Stan: A probabilistic programming language**. *Journal of statistical software*, 76(1).

