#ifndef FASTBC_SUBGRAPH_H
#define FASTBC_SUBGRAPH_H

#include "IGraph.h"
#include "ISubGraph.h"

#include <algorithm>
#include <memory>
#include <vector>

namespace fastbc {

	template<typename V, typename W>
	class SubGraph : public ISubGraph<V, W>
	{
	public:
		/**
		 *	@brief Initialize a sub-graph with given vertices
		 * 
		 *	@details Initialization will search for border vertices in O(m*log(n)) time
		 *			 where n is number of vertices and m number of edges
		 * 
		 *	@param subGraphVertices Vertices of sub-graph to consider
		 *	@param referenceGraph Full graph where the sub-graph is computed
		 */
		SubGraph(
			const std::vector<V>& subGraphVertices, 
			std::shared_ptr<const IGraph<V, W>> referenceGraph);

		W edge(V src, V dest) const override;

		const std::map<V, W>& forwardStar(V src) const override;

		const std::map<V, W>& backwardStar(V dest) const override;

		const std::vector<V>& verticesList() const override;

		V vertices() const override;

		V edges() const override;

		std::vector<V> borderVertices() const override;

		bool isBorder(V vertex) const override;

		std::shared_ptr<const IGraph<V, W>> referenceGraph() const override;

	private:
		const std::shared_ptr<const IGraph<V, W>> _referenceGraph;
		const std::vector<V> _vertices;
		V _edges;
		std::map<V, std::map<V, W>> _borderDestWeight;
		std::map<V, std::map<V, W>> _borderSrcWeight;
	};

}

template<typename V, typename W>
fastbc::SubGraph<V, W>::SubGraph(
	const std::vector<V>& subGraphVertices, 
	std::shared_ptr<const IGraph<V, W>> referenceGraph)
	: _referenceGraph(referenceGraph),
	_vertices(subGraphVertices),
	_edges(0)
{
	for (auto& v : _vertices)
	{
		// Check vertex forward star for edges terminating outside the graph
		const auto& fs = _referenceGraph->forwardStar(v);

		bool isBorder = false;
		std::vector<V> outEdges;

		for (auto& e : fs)
		{
			// When a vertex has an edge outside the sub-graph set it as border and store outgoing edge
			if (auto dest = std::find(_vertices.begin(), _vertices.end(), e.first); dest == _vertices.end())
			{
				isBorder = true;
				outEdges.push_back(e.first);
			}
		}

		// If vertex has been detected as border store a consistent forward star contained in sub-graph
		if (isBorder)
		{
			// Copy forward star from reference graph removing each unnecessary outgoing edge
			_borderDestWeight[v] = std::map<V, W>(fs);
			for (auto& out : outEdges)
			{
				_borderDestWeight[v].erase(out);
			}
		}

		// Update sub-graph edges counter
		_edges += fs.size() - outEdges.size();

		// Reset border variable for backward star check
		isBorder = false;
		outEdges.clear();

		// Check backward star for edges coming from outside the graph
		const auto& bs = _referenceGraph->backwardStar(v);
		for (auto& e : bs)
		{
			if (auto src = std::find(_vertices.begin(), _vertices.end(), e.first); src == _vertices.end())
			{
				isBorder = true;
				outEdges.push_back(e.first);
			}
		}

		if (isBorder)
		{
			_borderSrcWeight[v] = std::map<V, W>(bs);
			for (auto& out : outEdges)
			{
				_borderSrcWeight[v].erase(out);
			}
		}
	}
}

template<typename V, typename W>
W fastbc::SubGraph<V, W>::edge(V src, V dest) const
{
	const auto& fs = forwardStar(src);

	if (auto w = fs.find(dest); w != fs.end())
	{
		return w->second;
	}
	else
	{
		return 0;
	}
}

template<typename V, typename W>
const std::map<V, W>& fastbc::SubGraph<V, W>::forwardStar(V src) const
{
	if (auto border = _borderDestWeight.find(src); border != _borderDestWeight.end())
	{
		return border->second;
	}
	else
	{
		return _referenceGraph->forwardStar(src);
	}
}

template<typename V, typename W>
const std::map<V, W>& fastbc::SubGraph<V, W>::backwardStar(V dest) const
{
	if (auto border = _borderSrcWeight.find(dest); border != _borderSrcWeight.end())
	{
		return border->second;
	}
	else
	{
		return _referenceGraph->backwardStar(dest);
	}
}

template<typename V, typename W>
const std::vector<V>& fastbc::SubGraph<V, W>::verticesList() const
{
	return _vertices;
}

template<typename V, typename W>
V fastbc::SubGraph<V, W>::vertices() const
{
	return _vertices.size();
}

template<typename V, typename W>
V fastbc::SubGraph<V, W>::edges() const
{
	return _edges;
}

template<typename V, typename W>
std::vector<V> fastbc::SubGraph<V, W>::borderVertices() const
{
	std::vector<V> borders;

	for (auto& bdw : _borderDestWeight)
	{
		borders.push_back(bdw.first);
	}

	return borders;
}

template<typename V, typename W>
bool fastbc::SubGraph<V, W>::isBorder(V vertex) const
{
	return _borderDestWeight.find(vertex) != _borderDestWeight.end();
}

template<typename V, typename W>
std::shared_ptr<const fastbc::IGraph<V, W>> fastbc::SubGraph<V, W>::referenceGraph() const
{
	return _referenceGraph;
}

#endif