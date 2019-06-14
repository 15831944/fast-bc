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
make
```

### Usage
```
fbc [ options ] <edge_list_path>
```
```edge_list_path``` is a file containing the input graph represented as a list of edges. Each line represents an edge in the form ```<src> <dst> <weight>```. Each vertex is represented as an integer value (0 to #Vertices-1) and the graph is assumed to be directed with positive integer weights. 

For example, the following graph:

![](https://www.geeksforgeeks.org/wp-content/uploads/graph11.png)

can be represented with the following list of edges:
```
0 1 10
0 2 3
0 3 2
1 3 7
2 3 6
```

The output is a list of values where the value in position i is the betweennes centrality of the i-th vertex.

### Parameters

|Option   |Default value|Info|
|---|---|---|
|-s<br>--louvain-seeds||Louvain is an euristic algorithm. The output depends on the random order in which vertexes are examined. With this option you can pass a seed (int) to each louvain instance, to ensure the repeatability of results.|
|-e<br>--louvain-instances|4|To get better results, for each iteration of the Louvain algorithm the communities are calculated multiple times in parallel. In each parallel instance a different order for vertices examination is considered. The result with better modularity is then kept for the next iteraton. This parameter specify how many parallel instances of the partition calculation must run at each iteration.|
|-p<br>--louvain-precision|0.01|Terminate the Louvain algorithm when the difference in modularity between consecutive iterations is less than ```louvain-precision```.|
|  <br>--exact| |Force exact betweenness computation
|-t<br>--threads|OMP_NUM_THREADS|Maximum number of threads used in parallel computation|
|-k<br>--kfrac||Specify the number of superclasses that the second level of clustering must create. If for example, inside Louvain community 0 there are 100 classes and kfrac=0.5, the second level of clustering (kmeans) will generate 50 superclasses. |
|-o<br>--output|bc.txt|The output file name.|
|-d<br>--debug|info|Logger level (trace\|debug\|info\|warning\|error\|critical\|off)|

## References

|Algorithm   |Paper   |
|---|---|
|fast-bc   |A.Furno, N. El-Faouzi, R. Sharma, E. Zimeo. *Fast Computation of Betweenness Centrality to Locate Vulnerabilities in Very Large Road Networks*. In 97th Annual Meeting of the Transportation Research Board, July 2017. https://people.licit-lyon.eu/furno/documents/furno_trb18.pdf|
|Louvain   |V.D. Blondel,J.L. Guillaume,R. Lambiotte,E. Lefebvre. *Fast unfolding of communities in large networks.*. JSTAT 2008: P10008. https://arxiv.org/pdf/0803.0476.pdf   |
|Brandes   |U. Brandes.*A faster algorithm for betweenness centrality*. Journal of Mathematical Sociol-16ogy, 25(163), 2001. https://kops.uni-konstanz.de/bitstream/handle/123456789/5739/algorithm.pdf   |
