# Fast Computation of Betweenness Centrality

An implementation of [A Clustered Approach for Fast Computation of Betweenness Centrality](https://people.licit-lyon.eu/furno/documents/furno_trb18.pdf). The algorithm was adapted to analyze directed-weighted networks.

## Requirements
* C++17
* [CMake](https://cmake.org/) (>= 3.12)
* OpenMP (>= 4.5)
* [spdlog](https://github.com/gabime/spdlog) (>= 1.3.1)

## Getting Started

### Compile

```
git clone https://github.com/Daddeee/fast-bc
cd fast-bc/build
cmake ..
```

### Run

```
Usage: fbc [ options ] <edge_list_path>:
  -s, --louvain-seeds arg               Seeds to be used by each parallel louvain execution
  -e, --louvain-instances arg (=4)      Number of louvain instances
  -p, --louvain-precision arg (=0.01)   Minimum precision value for louvain algorithm
  - t -- threads arg (=OMP_NUM_THREADS) Number of threads.
  -k, --kfrac arg                       Topological classes aggregation factor (0-1)
  -o, --output arg (=bc.txt)            Output file path
  -d, --debug arg (=info)               Logger level (trace|debug|info|warning|error|critical|off)
  ```
  
## References

|Algorithm   |Paper   |
|---|---|
|fast-bc   |A.Furno, N. El-Faouzi, R. Sharma, E. Zimeo. *Fast Computation of Betweenness Centrality to Locate Vulnerabilities in Very Large Road Networks*. In 97th Annual Meeting of the Transportation Research Board, July 2017. https://people.licit-lyon.eu/furno/documents/furno_trb18.pdf|
|Louvain   |V.D. Blondel,J.L. Guillaume,R. Lambiotte,E. Lefebvre. *Fast unfolding of communities in large networks.*. JSTAT 2008: P10008. https://arxiv.org/pdf/0803.0476.pdf   |
|Brandes   |U. Brandes.*A faster algorithm for betweenness centrality*. Journal of Mathematical Sociol-16ogy, 25(163), 2001. https://kops.uni-konstanz.de/bitstream/handle/123456789/5739/algorithm.pdf   |
