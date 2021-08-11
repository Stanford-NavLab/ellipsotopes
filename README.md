# Ellipsotopes
This repository implements the ellipsotope, a set representation that combines the advantages of ellipsoids and zonotopes for reachability analysis and fault detection. Check out the paper [here](https://arxiv.org/abs/2108.01750).

The key aspect of this repository is the @ellipsotope class, which is implemented similar to ellipsoids in the [Ellipsoidal Toolbox](https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.html) and zonotopes in [CORA](https://github.com/TUMcps/CORA). The class implements a variety of operations (e.g., affine transformations, Minkowski sums, and intersections).

The "examples" folder contains the numerical examples on the paper, and the "figures" folder contains scripts to generate the figures in the paper. The "utilities" folder includes lots of helper functions needed by the rest of the code.

## Requirements
 - MATLAB (R2018a or greater)
 - [simulator](https://github.com/skousik/simulator) (for geometry and plotting utilities)
 - [MPT3](https://www.mpt3.org/) 
 - [CORA](https://github.com/TUMcps/CORA)

## References
If you find this work useful, please cite our paper that is available [here](https://arxiv.org/abs/2108.01750).

Kousik, S., Dai, A. and Gao, G., 2021. Ellipsotopes: Combining Ellipsoids and Zonotopes for Reachability Analysis and Fault Detection. arXiv preprint arXiv:2108.01750.
