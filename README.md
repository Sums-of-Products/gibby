# Gibby

**Gibby** is an MCMC algorithm that samples directed acyclic graphs (DAGs) from the posterior distribution.  

## Algorithm

Gibby combines several moves to efficiently explore the DAG space:  

- **Fast Gibby** 
- **New Edge Reversal (REV)** 
- **Markov Blanket Resampling (MBR)** 

Additionally, a **pruning technique** can be integrated to discard low-scoring parent sets.

## Input

`sample.cpp` runs the algorithm, and contains instructions for adjusting parameters. The algorithm accepts either: 

- **Discrete data**, or  
- **Local-score files** using [GOBNILP](https://www.cs.york.ac.uk/aig/sw/gobnilp/) format (also for continuous or mixed data)

The folder `data` contains the datasets used in the experiments.

## Output

The algorithm outputs two files containing:  
- **Sampled DAG scores**  
- **Edge-probability matrix** (row = parent, columns = child)  

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
