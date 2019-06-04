#ifndef FASTBC_DIRECTED_WEIGHTED_GRAPH_H
#define FASTBC_DIRECTED_WEIGHTED_GRAPH_H

#include <IGraph.h>

#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace fastbc {

    template<typename V, typename W>
    class DirectedWeightedGraph : public IGraph<V, W>
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

        W edge(V src, V dest) const override;

        const std::map<V, W>& forwardStar(V src) const override;

		const std::map<V, W>& backwardStar(V dest) const override;

		const std::set<V>& vertices() const override;

        V edges() const override;
        
    private:
        V _edges;
		std::set<V> _vertices;
        std::vector<std::map<V, W>> _srcDestWeight;
		std::vector<std::map<V, W>> _destSrcWeight;
    };   

}

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

        if(_srcDestWeight.size() <= src)
        {
            _srcDestWeight.resize(src + 1);
        }

		if (_destSrcWeight.size() <= dest)
		{
			_destSrcWeight.resize(dest + 1);
		}

        _srcDestWeight[src][dest] = weight;
		_destSrcWeight[dest][src] = weight;
        _edges++;
    }

	// Ensure both forward and backward star containers share same size
	_destSrcWeight.resize(_srcDestWeight.size());

	// Initialize vertices list
	for (V v = 0; v < _srcDestWeight.size(); v++)
	{
		_vertices.insert(v);
	}
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
const std::set<V>& fastbc::DirectedWeightedGraph<V, W>::vertices() const
{
	return _vertices;
}

template<typename V, typename W>
V fastbc::DirectedWeightedGraph<V, W>::edges() const
{
    return _edges;
}

#endif