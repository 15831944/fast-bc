#ifndef FASTBC_H
#define FASTBC_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>

using namespace boost;

typedef adjacency_list<vecS, listS, directedS,
		property<vertex_index_t, int>,
		property<edge_capacity_t, double> > Graph;

typedef property_map< Graph, vertex_index_t>::type VertexIndexMap;

typedef std::vector<std::vector<int>> Partition;

#endif