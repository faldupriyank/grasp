# Domain-Specialized Cache Management for Graph Analytics
Source code for the cache management technique, GRASP, published in [Faldu et al., HPCA'20] and the evaluated benchmarks.

This repo contains two parts: (1) Benchmarks and (2) Simulator.

## Benchmarks
This repo implements a simple data structure optimization to improve cache efficiency for three Ligra applications -- PageRank, PageRankDelta and BellmanFord.
To make this repo standalone, it retains a large chunk of the code from the original Ligra framework along with changes from Balaji et al. [3] and Faldu et al. [2]. For detailed information on the Ligra framework, please consult their original [repository](https://github.com/jshun/ligra). To know how to apply vertex reordering to improve cache efficiency, please consult the repository for [DBG](https://github.com/faldupriyank/dbg).

Specifically, we merge two arrays to improve spatial locality for the above mentioned applications. Files with Opt suffix contains our optimized version of the applications.

## Simulator
This repo contains a simple trace-based simulation infrastructure contianing four cache management techniques -- LRU, Belady, PIN and Grasp. The simulator expects a file containing a header and a trace of L2 miss addresses as specified in read_file() of cache.h. Simulator is very simple, and assumes that all addresses are memory read. Simulator provides a first order effect of various policies on cache efficiency. However, it may not provide accurate results for write-dominated applications. 

Finally, note that in the paper [1], we used the Sniper simulator. We included the Sniper APIs used to instrument applications to pass the region bounds to the cache microarchitecture in the Ligra source code. Look for code under the macro _SNIPER_.

# How to build

## Benchmarks
```
export DBG_ROOT='directory where this repo is cloned'
cd ${DBG_ROOT}/apps
make clean; make cleansrc; make -j
```
For more details on how to run, see [DBG](https://github.com/faldupriyank/dbg).

## Simulator
```
export DBG_ROOT='directory where this repo is cloned'
cd ${DBG_ROOT}/trace-based-simulators
make clean; make POLICY=grasp;
```

# Sample Results
Miss-rate for various cache management techniques for five graph applications processing the web-Google dataset using the traces provided in the ${DBG_ROOT}/datasets directory. Datasets are reorderd using DBG. Simulation assumes 16-way associative 1MB cache.

| Miss-rate(%) |  LRU |  PIN | GRASP | Belady |
|--------------|:----:|:----:|:-----:|:------:|
| BC           | 60.2 | 62.4 |  36.5 |  33.3  |
| SSSP         | 92.1 | 91.3 |  85.7 |  61.3  |
| PR           | 87.9 | 75.7 |  64.7 |  56.5  |
| PRD          | 87.9 | 75.7 |  64.7 |  56.5  |
| Radii        | 84.2 | 72.4 |  62.7 |  51.0  |

# License and copyright of the code used from external repositories
The repo also contains code from multiple other repositories (e.g., [Ligra](https://github.com/jshun/ligra), [Graph-Reordering-IISWC18](https://github.com/CMUAbstract/Graph-Reordering-IISWC18), [GAP](https://github.com/sbeamer/gapbs), [SNIPER](http://snipersim.org), [DBG](https://github.com/faldupriyank/dbg)) and the original copyright and license constraints apply to their code.

# References
**Please cite the following if you use the source code from this repository in your research.**
```
@inproceedings{DBGFalduHPCA20,  
  author={Priyank Faldu and Jeff Diamond and Boris Grot},  
  booktitle={International Symposium on High-Performance Computer Architecture (HPCA)},  
  title="{Domain-Specialized Cache Management for Graph Analytics}",  
  year={2020},  
  month=feb,  
}

@inproceedings{DBGFalduIISWC19,  
  author={Priyank Faldu and Jeff Diamond and Boris Grot},  
  booktitle={International Symposium on Workload Characterization (IISWC)},  
  title="{A Closer Look at Lightweight Graph Reordering}",  
  year={2019},  
  month=nov,  
}
```
[1] P. Faldu and J. Diamond and B. Grot, "Domain-Specialized Cache Management for Graph Analytics", in *Proceedings of the International Symposium on High-Performance Computer Architecture (HPCA)*, February 2020.

[2] P. Faldu and J. Diamond and B. Grot, "A Closer Look at Lightweight Graph Reordering," in *Proceedings of the International Symposium on Workload Characterization (IISWC)*, November 2019.

[3] V. Balaji and B. Lucia, "When is Graph Reordering an Optimization? Studying the Effect of Lightweight Graph Reordering Across Applications and Input Graphs," in *Proceedings of the International Symposium on Workload Characterization (IISWC)*, September 2018.
