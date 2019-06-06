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
	// Compute graph partition using Louvain communities detection algorithm
	auto communities = _le->evaluateGraph(graph);

	// Global betweenness centrality storage
	std::valarray<W> globalBC((W)0, grpah->vertices().size());

	// Vertices topological information about their own cluster border vertices
	std::vector<std::shared_ptr<VertexInfo<V, W>>> verticesInfo(graph->vertices().size(), nullptr);

	// Pivot vertices
	std::vector<V> pivots;

	// Pivot class cardinality for each vertex
	std::valarray<V> verticesClassCardinality(graph->vertices().size(), 1);

	// For each detected community extract related sub-graph, evaluate it for internal BC
	// and perform topological analysis to get pivots and vertices class cardinality
	for (int i = 0; i < communities.size(); i++)
	{
		std::shared_ptr<ISubGraph<V, W>> cluster = std::make_shared<SubGraph<V, W>>(communities[i]->all(), graph);

		_ce.evaluateCluster(globalBC, verticesInfo, cluster);

		std::vector<V> clusterPivots = 
			_ps.selectPivots(globalBC, verticesInfo, verticesClassCardinality, cluster->vertices());

		pivots.insert(pivots.end(), clusterPivots.begin(), clusterPivots.end());
	}

	std::valarray<W> intraClusterBC(globalBC);

	// Compute pivot contribution
	for (auto& pivot : pivots)
	{
		std::valarray<W> pivotDependency = _ssb.singleSourceBrandes(pivot, graph);

		globalBC += (pivotDependency - intraClusterBC) * verticesClassCardinality;
	}

	return globalBC;
}

#endif