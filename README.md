# **Gibby**

**Gibby** is a highly optimized **Markov Chain Monte Carlo (MCMC)** algorithm for sampling **directed acyclic graphs (DAGs)** from their **posterior distribution**.  

## 🧩 Algorithm Overview

Gibby efficiently explores the vast space of DAGs using several complementary move types:

- **Fast Basic Moves (FBM):** Add, remove, or reverse a single edge.  
- **New Edge Reversal (REV):** Efficient edge reversal based on  
  [Castelo & Kočka (2009)](https://link.springer.com/article/10.1007/s10994-008-5057-7).  
- **Markov Blanket Resampling (MBR):** Large-scale structural resampling of Markov blankets,  
  following [Su & Borsuk (2016)](https://jmlr.org/papers/v17/su16a.html).  

A **score-pruning technique** can optionally be enabled to discard low-probability parent sets during local score computation, improving scalability and memory efficiency.

## ⚙️ Compilation and Execution

Gibby is implemented in **C++17** and can be compiled with **g++** as follows:

```bash
g++ -std=c++17 -march=native -O3 -o gibby gibby.cpp
```

The mandatory parameters are the data file (the datasets used in the experiments are  in the data folder) and the number of main iterations of the MCMC algorithm. For example, run

```bash
./gibby  data/asia1k.dat -iter 1000
```
## 🔧 Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-d maxind` | Maximum in-degree per node | `-1` (unrestricted) |
| `-e ess` | Equivalent sample size (ESS) parameter for BDeu score | `1.0` |
| `-p prior` | Structure prior: `0` uniform, `1` fair, `2` fair+, `-w` edge(w) | `1` |
| `-K max` | Maximum number of candidate parents per node | `0` (unrestricted) |
| `-P mode` | Score pruning mode: `0` none, `1` top-down, `2` bottom-up; if preceded by digit `k>0`, only up to `min{k, maxind}` parents | `0` |
| `-a accuracy` | Number of significant bits in approximations | `15` |
| `-R seed` | Seed for random number generator | random |
| `-O file` | Output file name for local scores (in `.jkl` format) | none |
| `-M mem` | Amount of RAM available (in GiB) | `16` |
| `-burn in` | Number of burn-in iterations | `iter / 10` |
| `-FBM` | Number of Fast Basic Move steps per main iteration | `10000` |
| `-REV` | Number of REV moves per main iteration | `200` |
| `-MBR` | Number of MBR moves per main iteration | `200` |

## Output

The algorithm outputs two files containing:  
- **Sampled DAGs scores**  
- **Edge-probability matrix** (row = parent, column = child)

It is also possible to output another file containing the computed local scores (if the parameter -O is given).
