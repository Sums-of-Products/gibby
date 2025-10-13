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

Optional parameters are 

        -d maxind		                The maximum indegree parameter (default -1, i.e., infinity)
        -e ess		                    The ESS parameter of BDeu (default 1.0)
        -p prior		                The structure prior: 0 uniform, 1 fair, 2 fair+, -w edge(w) (default 1)
        -K max		                    Maximum number of parents per node (default 0 = unrestricted)
        -P mode		                    Score pruning: 0 no pruning, 1 top-down, 2 bottom-up (default 0); 
                                        if preceded by digit k > 0, then only up to min{k, maxind} parents
        -a accuracy		                The number of significant bits in approximations (default 15)
        -R seed		                    Seed to the random number generator (default random)
        -O filen		                Output file name for local scores in the jkl format 
        -M mem		                    Amount of RAM available in GiB (default 16)   
        
        -burn in                        number of main burn-in iterations (default iter/10)
        -FBM                            number of fast basic moves per main iteration (default 10000)
        -REV                            number of REV moves per main iteration (default 200)
        -MBR                            number of MBR moves per main iteration (default 200)
        
## Input

`sample.cpp` runs the algorithm, and contains instructions for adjusting parameters. The algorithm accepts either: 

- **Discrete data**, or  
- **Local-score files** using [GOBNILP](https://www.cs.york.ac.uk/aig/sw/gobnilp/) output format (also for continuous or mixed data)

The folder `data` contains the datasets used in the experiments.

## Output

The algorithm outputs two files containing:  
- **Sampled DAGs scores**  
- **Edge-probability matrix** (row = parent, column = child)

It is also possible to output another file containing the computed local scores (if the parameter -O is given).
