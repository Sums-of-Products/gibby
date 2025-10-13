# Gibby

**Gibby** is an highly optimized MCMC algorithm that samples directed acyclic graphs (DAGs) from the posterior distribution.  

## Algorithm

Gibby combines several moves to efficiently explore the space of DAGs:  

- **Fast basic moves** (add, remove, or reverse a single edge)
- [New Edge Reversal](https://link.springer.com/article/10.1007/s10994-008-5057-7) **(REV)** 
- [Markov Blanket Resampling](https://jmlr.org/papers/v17/su16a.html) **(MBR)** 

Additionally, a **pruning technique** can be integrated to discard low-scoring parent sets during the computation of local scores.

## Compilation, selection of parameters, and execution

This project is implemented in **C++**. To compile `sample.cpp` using **g++**, run

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

It is also possible to output another file containing the computed local scores. 

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
