# Gibby

Gibby is a MCMC algorithm that samples DAGs from the posterior distribution. Inputs can be either discrete data or local-score files in GOBNILP format. The algorithm combines fast Gibby, New Edge Reversal (REV), and Markov Blanket Resampling (MBR) moves, with pruning to eliminate low-scoring parent sets. 

`sample.cpp` runs the algorithm, outputs sampled DAG scores and a matrix of edge probabilities (parent = row, child = column), and contains instructions for adjusting parameters. The folder `data` contains the datasets used in the experiments.

## Language

This project is implemented in **C++**.

## Compilation

To compile the program using **g++**, run:

```bash
g++ -std=c++17 -march=native -O3 -o sample sample.cpp```

Once compiled, run using command

```bash
./sample
