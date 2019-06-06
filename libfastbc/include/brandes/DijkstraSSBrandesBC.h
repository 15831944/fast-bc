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
			std::valarray<W> singleSourceBrandes(
				V source,
				std::shared_ptr<const IGraph<V, W>> graph) override;

		private:

			template<typename N, typename E>
			struct vertex_backtrack_info_t
			{
				E sigma = 0;
				std::list<N> spPred;
			};

			template<typename N, typename E>
			struct backtrack_info_t
			{
				std::stack<N> visitStack;
				std::vector<vertex_backtrack_info_t<N, E>> spBacktrack;
			};

			backtrack_info_t<V, W> _dijkstra_SSSP(
				V src,
				std::shared_ptr<const IGraph<V, W>> graph);
		};

	}
}

template<typename V, typename W>
std::valarray<W> fastbc::brandes::DijkstraSSBrandesBC<V, W>::singleSourceBrandes(
	V source,
	std::shared_ptr<const IGraph<V, W>> graph)
{
	// Compute shortest path storing border information 
	backtrack_info_t<V, W> bi = _dijkstra_SSSP(source, cluster);
	auto& visitStack = bi.visitStack;
	auto& backtrackInfo = bi.spBacktrack;

	// Partial vertices dependency map
	std::vector<W> delta(cluster->vertices(), 0);

	std::valarray<W> ssBC(graph->vertices().size());

	// Backward visit of each vertex from dijkstra iteration 
	while (!visitStack.empty())
	{
		V w = visitStack.top();
		visitStack.pop();

		// Compute each vertex dependency for current src
		for (auto& v : backtrackInfo[w].spPred)
		{
			W c = backtrackInfo[v].sigma / backtrackInfo[w].sigma * (1.0 + delta[w]);

			delta[v] += c;
		}

		if (w != src)
		{
			ssBC[w] += delta[w];
		}
	}

	return ssBC;
}

template<typename V, typename W>
fastbc::brandes::DijkstraSSBrandesBC<V, W>::backtrack_info_t<V, W>
fastbc::brandes::DijkstraSSBrandesBC<V, W>::_dijkstra_SSSP(
	V src,
	std::shared_ptr<const IGraph<V, W>> graph)
{
	// Output information data structure
	backtrack_info_t<V, W> backtrackInfo;
	auto& visitStack = backtrackInfo.visitStack;
	auto& vertexBInfo = backtrackInfo.spBacktrack;
	vertexBInfo.resize(graph->vertices().size());

	// Map holding distances from the source.
	std::vector<W> dist(graph->vertices(), std::numeric_limits<W>::max());

	// Queue used for the Dijkstra's algorithm. Ordered by nearest vertex
	std::set<std::pair<W, V>> visitQueue;

	// Init src information
	vertexBInfo[src].sigma = 1;
	dist[src] = 0;
	visitQueue.insert(std::make_pair(dist[src], src));

	// While there are still elements in the queue.
	while (!visitQueue.empty())
	{
		// Pop the first.
		auto vPair = visitQueue.begin();
		V v = vPair->second;
		W srcDist = vPair->first;
		visitQueue.erase(vPair);

		// Push vertex to visited stack
		visitStack.push(v);

		// Check the neighbors w of v.
		for (const auto& it : graph->forwardStar(v))
		{
			V w = it.first;
			W newDist = srcDist + it.second;

			// Node w found for the first time or the new distance is shorter?
			if (newDist < dist[w])
			{
				visitQueue.erase(std::make_pair(dist[w], w));
				visitQueue.insert(std::make_pair(newDist, w));
				dist[w] = newDist;
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