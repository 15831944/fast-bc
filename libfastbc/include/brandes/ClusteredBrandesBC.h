#ifndef FASTBC_BRANDES_CLUSTEREDBRANDESBC_H
#define FASTBC_BRANDES_CLUSTEREDBRANDESBC_H

#include "IBrandesBC.h"
#include "IClusterEvaluator.h"
#include "ISSBrandesBC.h"
#include "IPivotSelector.h"
#include "VertexInfo.h"
#include <louvain/ILouvainEvaluator.h>
#include <SubGraph.h>

#include <memory>
#include <spdlog/spdlog.h>
#include <valarray>
#include <vector>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class ClusteredBrandeBC : public IBrandesBC<V, W>
		{
		public:
			ClusteredBrandeBC(
				std::shared_ptr<louvain::ILouvainEvaluator<V, W>> le,
				std::shared_ptr<IClusterEvaluator<V, W>> ce,
				std::shared_ptr<ISSBrandesBC<V, W>> ssb,
				std::shared_ptr<IPivotSelector<V, W>> ps);

			std::valarray<W> computeBC(const std::shared_ptr<const IGraph<V, W>> graph) override;

		private:
			std::shared_ptr<louvain::ILouvainEvaluator<V, W>> _le;
			std::shared_ptr<IClusterEvaluator<V, W>> _ce;
			std::shared_ptr<ISSBrandesBC<V, W>> _ssb;
			std::shared_ptr<IPivotSelector<V, W>> _ps;
		};

	}
}

template<typename V, typename W>
fastbc::brandes::ClusteredBrandeBC<V, W>::ClusteredBrandeBC(
	std::shared_ptr<fastbc::louvain::ILouvainEvaluator<V, W>> le,
	std::shared_ptr<fastbc::brandes::IClusterEvaluator<V, W>> ce,
	std::shared_ptr<fastbc::brandes::ISSBrandesBC<V, W>> ssb,
	std::shared_ptr<fastbc::brandes::IPivotSelector<V, W>> ps)
	: _le(le), _ce(ce), _ssb(ssb), _ps(ps)
{
}

template<typename V, typename W>
std::valarray<W> fastbc::brandes::ClusteredBrandeBC<V, W>::computeBC(
	const std::shared_ptr<const fastbc::IGraph<V, W>> graph)
{
	// Global betweenness centrality storage
	std::valarray<W> globalBC((W)0, graph->vertices().size());

	// Vertices topological information about their own cluster border vertices
	std::vector<std::shared_ptr<VertexInfo<V, W>>> verticesInfo(graph->vertices().size(), nullptr);

	// Computed subgraph and border vertices from each vertices community
	std::vector<std::shared_ptr<ISubGraph<V, W>>> cluster;

	// Pivot vertices for each cluster
	std::vector<std::vector<V>> pivotsCluster;

	// Pivot class cardinality for each vertex
	std::valarray<W> verticesClassCardinality((W)1, graph->vertices().size());

	// Compute graph partition using Louvain communities detection algorithm
	SPDLOG_INFO("Computing clusters with Louvain algorithm");
	std::vector<std::shared_ptr<louvain::ICommunity<V, W>>> communities = 
		_le->evaluateGraph(std::static_pointer_cast<const IDegreeGraph<V, W>>(graph));

	SPDLOG_INFO("Graph partitioned in {} clusters", communities.size());
	cluster.resize(communities.size());
	pivotsCluster.resize(communities.size());

	// For each detected community compute related sub-graph, evaluate it for internal BC
	// and perform topological analysis to get pivots and vertices class cardinality
	for (int i = 0; i < cluster.size(); i++)
	{
		cluster[i] = std::make_shared<SubGraph<V, W>>(communities[i]->all(), graph);

		SPDLOG_INFO("Evaluating BC on cluster {}: {} vertices, {} edges", 
			i, cluster[i]->vertices().size(), cluster[i]->edges());
		_ce->evaluateCluster(globalBC, verticesInfo, cluster[i]);

		pivotsCluster[i] = _ps->selectPivots(
			globalBC, verticesInfo, verticesClassCardinality, 
			cluster[i]->vertices(), cluster[i]->borders());

		SPDLOG_INFO("Selected {} vertices as pivots in cluster {}", pivotsCluster[i].size(), i);
	}

	std::valarray<W> intraClusterBC(globalBC);

	
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_INFO
	size_t pivotCount = 0;
	for (const auto& pc : pivotsCluster)
	{
		pivotCount += pc.size();
	}
	SPDLOG_INFO("Computing global BC from {} pivots", pivotCount);
#endif

	// Compute global dependecy pivot contribution
	for (size_t c = 0; c < cluster.size(); ++c)
	{
		for (size_t p = 0; p < pivotsCluster[c].size(); ++p)
		{
			std::valarray<W> pivotDependency = _ssb->singleSourceBrandes(pivotsCluster[c][p], graph);

			// Sum pivot dependecy to all vertices
			globalBC += pivotDependency * verticesClassCardinality[pivotsCluster[c][p]];

			// Subtract duplicate dependency from current pivot's cluster vertices
			for (const auto& v : cluster[c]->vertices())
			{
				globalBC[v] -= intraClusterBC[v] * verticesClassCardinality[pivotsCluster[c][p]];
			}
		}
	}

	return globalBC;
}

#endif