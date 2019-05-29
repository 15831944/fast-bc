#ifndef FASTBC_DIRECTED_WEIGHTED_GRAPH_H
#define FASTBC_DIRECTED_WEIGHTED_GRAPH_H

#include <IGraph.h>

#include <iostream>
#include <map>
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
		 *	@param inputTextGraph Graph edges input stream
		 */
        DirectedWeightedGraph(std::istream& inputTextGraph);

        W edge(V src, V dest) const;

        const std::map<V, W>& forwardStar(V src) const;

        V nodes() const ;

        V edges() const;
        
    private:
        V _edges;
        std::vector<std::map<V, W>> _srcDestWeight;
    };   

}

template<typename V, typename W>
fastbc::DirectedWeightedGraph<V, W>::DirectedWeightedGraph(std::istream& inputTextGraph)
    : _edges(0), _srcDestWeight(10)
{
    while (!inputTextGraph.eof())
    {
        V src, dest;
		W weight;
        inputTextGraph >> src >> dest >> weight;

        if(_srcDestWeight.size() <= src)
        {
            _srcDestWeight.resize(src + 1);
        }

        _srcDestWeight[src][dest] = weight;
        _edges++;
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
V fastbc::DirectedWeightedGraph<V, W>::nodes() const
{
    return _srcDestWeight.size();
}

template<typename V, typename W>
V fastbc::DirectedWeightedGraph<V, W>::edges() const
{
    return _edges;
}

#endif