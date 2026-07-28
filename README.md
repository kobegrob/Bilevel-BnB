# Bilevel Branch-and-Bound Solver

This repository contains an implementation of the combinatorial branch-and-bound algorithm introduced in *(insert archive reference here)* for solving bilevel mixed-integer packing problems with interdiction constraints.

> **Note:** The current version has only been tested on **Linux**.

---

## Building

To compile the code, run the following command in your terminal:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -o bilevel Bilevel.cpp
```
After creating the executable, it can be called using: 

```bash
./bilevel
```
You will be promted to give a file name (.knp file) without the extention.

Alternatively, it is possible to run a folder of files and save the results in a .csv file:
```bash
./run_bilevel_batch.sh foldername filename.csv
```
## Input file
The code expects .knp files which are structured as follows:
The first line contains n  C^U  C^L (which represent the number of items and the upper- and lower-level capacities). 
The next n lines contain w  f  c  l  v (which are the parameters explained in the paper for each of the items).

## Overview of the results
This repository contains all experimental results presented in the accompanying paper. The proposed combinatorial branch-and-bound (BnB) approach is benchmarked against the following state-of-the-art bilevel optimization solvers:

- The General-purpose bilevel solver by Fischetti et al. (2017) (https://msinnl.github.io/pages/bilevel.html) (Results can be found in 'Results General Solver' folder).
  
- MiBs solver by Tahernejad et al. (2020) (https://github.com/coin-or/MibS) (Results can be found in 'Results MiBs' folder).

All instance files for used for the branch-and-bound approach can be found in the folder 'Instance files'. The'Results BnB' folder contains all computational results obtained with the proposed branch-and-bound method and reported in the paper. 

The most important results are: 

- 

### References
- Fischetti, M., Ljubić, I., Monaci, M., & Sinnl, M. (2017). A New General-Purpose Algorithm for Mixed-Integer Bilevel Linear Programs. Operations Research, 65(6), 1615–1637. https://doi.org/10.1287/opre.2017.1650
- Tahernejad, S., Ralphs, T. K., & DeNegre, S. T. (2020). A branch-and-cut algorithm for mixed integer bilevel linear optimization problems and its implementation. Mathematical Programming Computation, 12(4), 529–568. https://doi.org/10.1007/s12532-020-00183-6



