#ifndef FASTBC_BRANDES_DIJKSTRACLUSTEREVALUATOR_H
#define FASTBC_BRANDES_DIJKSTRACLUSTEREVALUATOR_H

#include "IClusterEvaluator.h"

#include <list>
#include <memory>
#include <set>
#include <stack>
#include <tuple>
#include <vector>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class DijkstraClusterEvaluator : public IClusterEvaluator<V, W>
		{
		public:
			DijkstraClusterEvaluator();
			~DijkstraClusterEvaluator();

			void evaluateCluster(
				std::valarray<W>& clusterBC,
				std::vector<std::shared_ptr<VertexInfo<V, W>>>& globalVI,
				std::shared_ptr<const ISubGraph<V, W>> cluster) override;

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
				std::map<N, vertex_backtrack_info_t<N, E>> spBacktrack;
			};

			backtrack_info_t<V, W> dijkstra_SSSP(
				std::vector<std::shared_ptr<VertexInfo<V, W>>>& globalVI,
				V src,
				std::shared_ptr<const ISubGraph<V, W>> graph);

		};

	}
}

template<typename V, typename W>
fastbc::brandes::DijkstraClusterEvaluator<V, W>::DijkstraClusterEvaluator()
{
}

template<typename V, typename W>
fastbc::brandes::DijkstraClusterEvaluator<V, W>::~DijkstraClusterEvaluator()
{
}

template<typename V, typename W>
void fastbc::brandes::DijkstraClusterEvaluator<V, W>::evaluateCluster(
	std::valarray<W>& clusterBC,
	std::vector<std::shared_ptr<VertexInfo<V, W>>>& globalVI,
	std::shared_ptr<const ISubGraph<V, W>> cluster)
{
	// Partial vertices dependency map
	std::map<V, W> delta;
	for (const auto& v : cluster->vertices()) { delta[v] = 0; }

	// Compute SP from each cluster vertex
	for (const auto& src : cluster->vertices())
	{
		// Reset partial dependency structure before starting
		for (auto& vw : delta) { vw.second = 0; }

		// Compute shortest path storing border information 
		backtrack_info_t<V, W> bi = dijkstra_SSSP(globalVI, src, cluster);
		auto& visitStack = bi.visitStack;
		auto& backtrackInfo = bi.spBacktrack;

		// Backward visit of each vertex from dijkstra iteration 
		while (!bi.visitStack.empty())
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
				clusterBC[w] += delta[w];
			}
		}
	}
}

template<typename V, typename W>
fastbc::brandes::DijkstraClusterEvaluator<V, W>::backtrack_info_t<V, W>
fastbc::brandes::DijkstraClusterEvaluator<V, W>::dijkstra_SSSP(
	std::vector<std::shared_ptr<VertexInfo<V, W>>>& globalVI,
	V src,
	std::shared_ptr<const ISubGraph<V, W>> graph)
{
	// Output information data structure
	backtrack_info_t<V, W> backtrackInfo;
	auto& visitStack = backtrackInfo.visitStack;
	auto& vertexBInfo = backtrackInfo.spBacktrack;

	// Map holding distances from the source.
	std::map<V, W> dist;
	for (const auto& v : graph->vertices()) { dist[v] = std::numeric_limits<W>::max(); }

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

	// Annotate shortest path length and count information from current src to border vertices
	const auto& borders = graph->borders();
	V storeIndex = 0;
	globalVI[src] = std::make_shared<VertexInfo<V, W>>(borders.size());
	for (const auto& b : borders)
	{
		globalVI[src]->setBorderSPLength(storeIndex, dist[b]);
		globalVI[src]->setBorderSPCount(storeIndex, vertexBInfo[b].sigma);
		storeIndex++;
	}

	return backtrackInfo;
}

#endif