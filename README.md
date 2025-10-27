# **Gibby**

**Gibby** is a highly optimized **Markov Chain Monte Carlo (MCMC)** algorithm for sampling **directed acyclic graphs (DAGs)** from their **posterior distribution**.  

## Algorithm Overview

Gibby efficiently explores the vast space of DAGs using several move types:

- **Gibby Fast Basic Moves (FBM):** Adds, removes, or reverses a single edge.  
- **[New Edge Reversal Move](https://link.springer.com/article/10.1007/s10994-008-5057-7) (REV):** Reverses edges while resampling entire parent sets.
- **[Markov Blanket Resampling Move](https://jmlr.org/papers/v17/su16a.html) (MBR):** Performs large-scale resampling of Markov blankets.  

A **score-pruning technique** can optionally be enabled to discard low-probability parent sets during local score computation, improving scalability and memory efficiency.

## Compilation and Execution

Gibby is implemented in **C++** and can be compiled with **g++** as follows:

```bash
g++ -std=c++17 -march=native -O3 -o gibby gibby.cpp
```
The only mandatory parameter is the data file (the datasets used in the experiments are located in the `data` folder). For example, run the algorithm using:

```bash
./gibby  data/asia1k.dat 
```
Currently, Gibby can generate scores only from discrete data, scored using the BDeu system. 

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-d` | Maximum in-degree per node | `-1` (unrestricted) |
| `-e` | Equivalent sample size (ESS) parameter for BDeu score | `1.0` |
| `-p` | Structure prior: `0` uniform, `1` fair, `2` fair+, `-w` edge(w) | `1` |
| `-K` | Number of candidate parents per node | `0` (unrestricted) |
| `-P` | Score pruning mode: `0` none, `1` top-down, `2` bottom-up; if preceded by digit `k>0`, only up to `min{k, maxind}` parents | `0` |
| `-a` | Number of significant bits in approximations | `15` |
| `-R` | Seed for random number generator | random |
| `-O` | Output file name for local scores (in `.jkl` format) | none |
| `-M` | Amount of RAM available (in GiB) | `16` |
| `-iter` | Number of main iterations | `10 000` |
| `-burnin` | Number of burn-in iterations | `iter / 10` |
| `-FBM` | Number of Fast Basic Moves  per iteration | `10 000` |
| `-REV` | Number of REV moves per iteration | `200` |
| `-MBR` | Number of MBR moves per iteration | `200` |

### Run using precomputed local scores

For continuous or mixed networks, you can use precomputed local scores instead of raw data. These scores should be provided using the .jkl format. In this case, run the algorithm using:

```bash
./gibby  -I score_file_path.jkl
```

## Output Files

After execution, Gibby produces the following outputs:

- **Sampled DAGs scores** — contains the sampled DAGs scores over iterations.
- **Edge probability matrix** — contains the posterior probability of the edges  (row = parent, column = child).
- **(Optional) Parent sets scores** — if `-O <filename>` is specified, Gibby also outputs the computed local scores in `.jkl` format.
