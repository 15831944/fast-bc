#ifndef FASTBC_BRANDES_DIJKSTRASSBRANDESBC_H
#define FASTBC_BRANDES_DIJKSTRASSBRANDESBC_H

#include "ISSBrandesBC.h"

#include <list>
#include <set>
#include <stack>
#include <vector>
#include <utility>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class DijkstraSSBrandesBC : public ISSBrandesBC<V, W>
		{
		public:
			std::vector<W> singleSourceBrandes(
				V source,
				std::shared_ptr<const IGraph<V, W>> graph) override;

		private:

			struct vertex_backtrack_info_t
			{
				W sigma = 0;
				std::list<V> spPred;
			};

			struct backtrack_info_t
			{
				std::stack<V> visitStack;
				std::vector<vertex_backtrack_info_t> spBacktrack;
			};

			backtrack_info_t _dijkstra_SSSP(
				V src,
				std::shared_ptr<const IGraph<V, W>> graph);
		};

	}
}

template<typename V, typename W>
std::vector<W> fastbc::brandes::DijkstraSSBrandesBC<V, W>::singleSourceBrandes(
	V source,
	std::shared_ptr<const IGraph<V, W>> graph)
{
	// Compute shortest path storing border information 
	struct backtrack_info_t bi = _dijkstra_SSSP(source, graph);
	auto& visitStack = bi.visitStack;
	auto& backtrackInfo = bi.spBacktrack;

	// Partial vertices dependency
	std::vector<W> delta(graph->vertices().size(), 0);

	std::vector<W> ssBC(graph->vertices().size(), (W)0);

	// Backward visit of each vertex from dijkstra iteration 
	while (!visitStack.empty())
	{
		V w = visitStack.top();
		visitStack.pop();

		// Compute each vertex dependency for current src
		for (const auto& v : backtrackInfo[w].spPred)
		{
			W c = backtrackInfo[v].sigma / backtrackInfo[w].sigma * (1.0 + delta[w]);

			delta[v] += c;
		}

		if (w != source)
		{
			ssBC[w] += delta[w];
		}
	}

	return ssBC;
}

template<typename V, typename W>
struct fastbc::brandes::DijkstraSSBrandesBC<V, W>::backtrack_info_t
fastbc::brandes::DijkstraSSBrandesBC<V, W>::_dijkstra_SSSP(
	V src,
	std::shared_ptr<const IGraph<V, W>> graph)
{
	// Output information data structure
	struct backtrack_info_t backtrackInfo;
	auto& visitStack = backtrackInfo.visitStack;
	auto& vertexBInfo = backtrackInfo.spBacktrack;
	vertexBInfo.resize(graph->vertices().size());

	// Map holding distances from the source.
	std::vector<W> dist(graph->vertices().size(), std::numeric_limits<W>::max());

	// Queue used for the Dijkstra's algorithm. Ordered by nearest vertex to src
	auto distCmp = [&dist](const V& lhs, const V& rhs) { 
		if(dist[lhs] == dist[rhs])
			return lhs < rhs;
		return dist[lhs] < dist[rhs]; 
	};
	std::set<V, decltype(distCmp)> visitQueue(distCmp);

	// Init src information
	vertexBInfo[src].sigma = 1;
	dist[src] = 0;
	visitQueue.insert(src);

	// While there are still elements in the queue.
	while (!visitQueue.empty())
	{
		// Pop the first
		V v = *visitQueue.begin();
		visitQueue.erase(visitQueue.begin());

		// Push vertex to visited stack
		visitStack.push(v);

		// Check the neighbors w of v.
		for (const auto& it : graph->forwardStar(v))
		{
			V w = it.first;
			W newDist = dist[v] + it.second;

			// Node w found for the first time or the new distance is shorter?
			if (newDist < dist[w])
			{
				visitQueue.erase(w);
				dist[w] = newDist;
				visitQueue.insert(w);
				vertexBInfo[w].spPred.clear();
				vertexBInfo[w].sigma = 0;
			}

			// Is the shortest path to w via u?
			if (newDist == dist[w])
			{
				vertexBInfo[w].spPred.push_back(v);
				vertexBInfo[w].sigma += vertexBInfo[v].sigma;
			}
		}
	}

	return backtrackInfo;
}

#endif