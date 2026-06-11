# Bilevel Branch-and-Bound Solver

This repository contains an implementation of the combinatorial branch-and-bound algorithm introduced in *(insert paper reference here)* for solving bilevel mixed-integer packing problems with interdiction constraints.

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
This repository also contains all results from the experiments presented in the accompanying paper.
