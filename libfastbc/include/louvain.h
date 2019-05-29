#ifndef LOUVAIN_H
#define LOUVAIN_H

#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>

#include <fastbc.h>

class Louvain
{
public:
	Louvain(Graph g) : g(g) {}
	
	Partition compute(std::vector<int> evaluation_order) 
	{
		bool done;
		Graph g1 = g;
		Partition p1 = singleton_partition(g1);

		do {
			p1 = move_nodes(g1, p1);
			int size = num_vertices(g1);
			done = (p1.size() == size);
			if(!done) {
				g1 = aggregate_graph(g1, p1);
				p1 = singleton_partition(g1);
			} 
		} while(!done);

		return p1;
	}

	std::vector<int> random_evaluation_order()
	{
		std::vector<int> v;
		int size = num_vertices(g);
		boost::push_back(v, boost::irange(0, size));
		boost::range::random_shuffle(v);
		return v;
	}

private:
	Graph g;

	Partition move_nodes(Graph g1, Partition p1) 
	{
		do {
			double h_old = modularity(p1);
		} while(true);
	}

	Graph aggregate_graph(Graph g1, Partition p1) 
	{

	}

	Partition singleton_partition(Graph g1)
	{
		Partition p;
		VertexIndexMap v_index = get(vertex_index, g);
		BGL_FORALL_VERTICES_T(v, g1, Graph) {
			std::vector<int> vs(1, v_index[v]);
	    	p.push_back(vs);
		}
		return p;
	}

	double modularity(Partition p) {
		return 1.;
	}
};

#endif //LOUVAIN_H