# Gibby

**Gibby** is an MCMC algorithm that samples directed acyclic graphs (DAGs) from the posterior distribution.  

## Algorithm

Gibby combines several moves to efficiently explore the DAG space:  

- **Fast Gibby** 
- [New Edge Reversal](https://link.springer.com/article/10.1007/s10994-008-5057-7) **(REV)** 
- [Markov Blanket Resampling](https://jmlr.org/papers/v17/su16a.html) **(MBR)** 

Additionally, a **pruning technique** can be integrated to discard low-scoring parent sets.

## Input

`sample.cpp` runs the algorithm, and contains instructions for adjusting parameters. The algorithm accepts either: 

- **Discrete data**, or  
- **Local-score files** using [GOBNILP](https://www.cs.york.ac.uk/aig/sw/gobnilp/) format (also for continuous or mixed data)

The folder `data` contains the datasets used in the experiments.

## Output

The algorithm outputs two files containing:  
- **Sampled DAG scores**  
- **Edge-probability matrix** (row = parent, column = child)

It is possible to output a file containing the computed local scores. 

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
