#ifndef FASTBC_BRANDES_EXACTBRANDESBC_H
#define FASTBC_BRANDES_EXACTBRANDESBC_H

#include "IBrandesBC.h"

#include <functional>
#include <list>
#include <memory>
#include <set>
#include <stack>
#include <vector>

namespace fastbc {
    namespace brandes {

        template<typename V, typename W>
        class ExactBrandesBC : public IBrandesBC<V, W>
        {
        public:
            std::vector<W> computeBC(const std::shared_ptr<const IGraph<V, W>> graph) override;


        private:

            struct vertex_backtrack_info_t
			{
				W sigma = 0;
				std::list<V> spPred;
			};

			struct backtrack_info_t
			{
				std::stack<V> visitStack;
				std::vector<struct vertex_backtrack_info_t> spBacktrack;
			};

			backtrack_info_t _dijkstra_SSSP(
				V src,
				std::shared_ptr<const IGraph<V, W>> graph);
        };

    }
}

template<typename V, typename W>
std::vector<W> fastbc::brandes::ExactBrandesBC<V, W>::computeBC(
    const std::shared_ptr<const IGraph<V, W>> graph)
{
    std::vector<W> globalBC(graph->vertices().size(), (W)0);
    W* _globalBC = globalBC.data();
	size_t _globalBCsize = globalBC.size();

	#pragma omp parallel
	{
		// Partial dependency vertices map
		std::vector<W> delta(graph->vertices().size(), (W)0);

		// Compute SP from each cluster vertex
		#pragma omp for reduction(+:_globalBC[:_globalBCsize])
		for (size_t srcIndex = 0; srcIndex < graph->vertices().size(); ++srcIndex)
		{
			const V& src = graph->vertices()[srcIndex];

			// Reset partial dependency structure before starting
			delta.assign(delta.size(), 0);

			// Compute shortest path storing border information 
			struct backtrack_info_t bi = _dijkstra_SSSP(src, graph);
			auto& visitStack = bi.visitStack;
			auto& backtrackInfo = bi.spBacktrack;

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
					_globalBC[w] += delta[w];
				}
			}
		}
	}

    return globalBC;
}

template<typename V, typename W>
struct fastbc::brandes::ExactBrandesBC<V, W>::backtrack_info_t
fastbc::brandes::ExactBrandesBC<V, W>::_dijkstra_SSSP(
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