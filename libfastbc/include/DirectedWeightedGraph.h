#ifndef FASTBC_DIRECTED_WEIGHTED_GRAPH_H
#define FASTBC_DIRECTED_WEIGHTED_GRAPH_H

#include <IDegreeGraph.h>

#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace fastbc {

    template<typename V, typename W>
    class DirectedWeightedGraph : public IDegreeGraph<V, W>
    {
    public:
		/**
		 *	@brief Initialize a directed weighted graph from input stream
		 *	
		 *	@details Given input stream should feed edges information like:
		 *			<src_index> <dest_index> <wdge_weight>
		 * 
		 *	@note Input graph is expected to have subsequent vertex indices starting from 0
		 *
		 *	@param inputTextGraph Graph edges input stream
		 */
        DirectedWeightedGraph(std::istream& inputTextGraph);

        DirectedWeightedGraph();

        W edge(V src, V dest) const override;

        const std::map<V, W>& forwardStar(V src) const override;

		const std::map<V, W>& backwardStar(V dest) const override;

		const std::vector<V>& vertices() const override;

        V edges() const override;

        void addEdge(V from, V to, W weight) override;

        void initVertices() override;

        W totalWeight() const override;

		W inWeightedDegree(V v) const override;

        W outWeightedDegree(V v) const override;
        
    private:
        V _edges;
        W _totalWeight;
		std::vector<V> _vertices;
		std::vector<W> _inWeightedDegrees;
		std::vector<W> _outWeightedDegrees;
        std::vector<std::map<V, W>> _srcDestWeight;
		std::vector<std::map<V, W>> _destSrcWeight;
    };   

}

template<typename V, typename W>
fastbc::DirectedWeightedGraph<V, W>::DirectedWeightedGraph()
	: _edges(0) {}

template<typename V, typename W>
fastbc::DirectedWeightedGraph<V, W>::DirectedWeightedGraph(std::istream& inputTextGraph)
    : _edges(0)
{
	// Read input stream and initialize forward and backward star for each vertex
    while (!inputTextGraph.eof())
    {
        V src, dest;
		W weight;
        inputTextGraph >> src >> dest >> weight;

		if (weight <= 0)
		{
			throw std::invalid_argument("Edge weight must be greater than zero");
		}

		// Stop when eof has been reached
		if (!inputTextGraph)
		{
			break;
		}

        addEdge(src, dest, weight);
    }

	// Ensure both forward and backward star containers share same size
	_destSrcWeight.resize(_srcDestWeight.size());

	initVertices();
}

template<typename V, typename W>
W fastbc::DirectedWeightedGraph<V, W>::edge(V src, V dest) const
{
    if(auto it = _srcDestWeight[src].find(dest); it != _srcDestWeight[src].end())
    {
        return it->second;
    }
    else
    {
        return 0;
    }
}

template<typename V, typename W>
const std::map<V, W>& fastbc::DirectedWeightedGraph<V, W>::forwardStar(V src) const
{
    return _srcDestWeight[src];
}

template<typename V, typename W>
const std::map<V, W>& fastbc::DirectedWeightedGraph<V, W>::backwardStar(V dest) const
{
	return _destSrcWeight[dest];
}

template<typename V, typename W>
const std::vector<V>& fastbc::DirectedWeightedGraph<V, W>::vertices() const
{
	return _vertices;
}

template<typename V, typename W>
V fastbc::DirectedWeightedGraph<V, W>::edges() const
{
    return _edges;
}

template<typename V, typename W>
void fastbc::DirectedWeightedGraph<V, W>::addEdge(V from, V to, W weight) 
{
	if(_srcDestWeight.size() <= from)
    {
        _srcDestWeight.resize(from + 1);
    }

	if (_destSrcWeight.size() <= to)
	{
		_destSrcWeight.resize(to + 1);
	}

	int s = (to > from) ? to : from;
	if (_inWeightedDegrees.size() <= s)
		_inWeightedDegrees.resize(s + 1, 0);
	if (_outWeightedDegrees.size() <= s)
		_outWeightedDegrees.resize(s + 1, 0);

	if(_srcDestWeight[from].count(to) > 0) {
		_srcDestWeight[from][to] += weight;
	} else {
		_srcDestWeight[from][to] = weight;
    	_edges++;
	}

	if(_destSrcWeight[to].count(from) > 0) {
		_destSrcWeight[to][from] += weight;
	} else {
		_destSrcWeight[to][from] = weight;
	}

    _totalWeight += weight;
	_inWeightedDegrees[to] += weight;
	_outWeightedDegrees[from] += weight;
}

template<typename V, typename W>
void fastbc::DirectedWeightedGraph<V, W>::initVertices() 
{
	// Initialize vertices list
	_vertices.resize(_srcDestWeight.size());
	#pragma omp simd
	for (size_t v = 0; v < _vertices.size(); v++)
	{
		_vertices[v] = v;
	}
}

template<typename V, typename W>
W fastbc::DirectedWeightedGraph<V, W>::totalWeight() const
{
	return _totalWeight;
}

template<typename V, typename W>
W fastbc::DirectedWeightedGraph<V, W>::inWeightedDegree(V v) const
{
	return _inWeightedDegrees[v];
}

template<typename V, typename W>
W fastbc::DirectedWeightedGraph<V, W>::outWeightedDegree(V v) const
{
	return _outWeightedDegrees[v];
}

#endif