# Gibby

**Gibby** is an MCMC algorithm that samples directed acyclic graphs (DAGs) from the posterior distribution.  

## Input

The algorithm accepts either:  
- **Discrete data**, or  
- **Local-score files** in **GOBNILP** format.  

## Algorithm

Gibby combines several moves to efficiently explore the DAG space:  

- **Fast Gibby** move
- **New Edge Reversal (REV)** move
- **Markov Blanket Resampling (MBR)** move

Additionally, a **pruning technique** is applied to discard low-scoring parent sets.

`sample.cpp` runs the algorithm, outputs sampled DAG scores and a matrix of edge probabilities (parent = row, child = column), and contains instructions for adjusting parameters. The folder `data` contains the datasets used in the experiments.

## Language

This project is implemented in **C++**.

## Compilation and execution

To compile the program using **g++**, run

```bash
g++ -std=c++17 -march=native -O3 -o sample sample.cpp
```

Once compiled, execute it using

```bash
./sample
```
