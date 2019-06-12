#ifndef FASTBC_SUBGRAPH_H
#define FASTBC_SUBGRAPH_H

#include "IGraph.h"
#include "ISubGraph.h"

#include <map>
#include <memory>
#include <set>
#include <spdlog/spdlog.h>
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

		const std::vector<V>& vertices() const override;

		V edges() const override;

		const std::set<V>& borders() const override;

		bool isBorder(V vertex) const override;

		std::shared_ptr<const IGraph<V, W>> referenceGraph() const override;

	private:
		const std::shared_ptr<const IGraph<V, W>> _referenceGraph;
		const std::vector<V> _vertices;
		V _edges;
		std::map<V, std::map<V, W>> _borderDestWeight;
		std::map<V, std::map<V, W>> _borderSrcWeight;
		std::set<V> _borderVertices;
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
	// Order sub-graph vertices to speed-up border vertices computation
	std::set orderedVertices(_vertices.begin(), _vertices.end());

	for (size_t vIndex = 0; vIndex < _vertices.size(); ++vIndex)
	{
		const V& v = _vertices[vIndex];

		// Check vertex forward star for edges terminating outside the graph
		const auto& fs = _referenceGraph->forwardStar(v);

		bool isBorder = false;
		V connections = 0;
		std::vector<V> outEdges;

		for (auto& e : fs)
		{
			// When a vertex has an edge outside the sub-graph set it as border and store outgoing edge
			if (auto dest = orderedVertices.find(e.first); dest == orderedVertices.end())
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
			connections += _borderDestWeight[v].size();

			_borderVertices.insert(v);
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
			if (auto src = orderedVertices.find(e.first); src == orderedVertices.end())
			{
				isBorder = true;
				outEdges.push_back(e.first);
			}
		}

		// If vertex has been detected as border store a consistent backward star contained in sub-graph
		if (isBorder)
		{
			// Copy backward star from reference graph removing each unnecessary external edge
			_borderSrcWeight[v] = std::map<V, W>(bs);
			for (auto& out : outEdges)
			{
				_borderSrcWeight[v].erase(out);
			}
			connections += _borderSrcWeight[v].size();

			_borderVertices.insert(v);
		}

		// If a vertex runs out of edges, the sub-graph is not consistent
		if (isBorder && !connections && !(_vertices.size() == 1))
		{
			SPDLOG_TRACE("Vertex {} is unconnected in its cluster", v);

#ifdef FASTBC_SUBGRAPH_CONNECTED_ONLY
			SPDLOG_CRITICAL("Vertex {} is unconnected in its cluster", v);
			throw std::invalid_argument("Given subgraph has unconnected vertices");
#endif
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
const std::vector<V>& fastbc::SubGraph<V, W>::vertices() const
{
	return _vertices;
}

template<typename V, typename W>
V fastbc::SubGraph<V, W>::edges() const
{
	return _edges;
}

template<typename V, typename W>
const std::set<V>& fastbc::SubGraph<V, W>::borders() const
{
	return _borderVertices;
}

template<typename V, typename W>
bool fastbc::SubGraph<V, W>::isBorder(V vertex) const
{
	return _borderVertices.find(vertex) != _borderVertices.end();
}

template<typename V, typename W>
std::shared_ptr<const fastbc::IGraph<V, W>> fastbc::SubGraph<V, W>::referenceGraph() const
{
	return _referenceGraph;
}

#endif