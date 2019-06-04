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
	std::valarray<W> globalBC((W)0, grpah->vertices());

	// Communities vertices information
	std::vector<std::shared_ptr<VertexInfo<V, W>>> commVertexInfos(graph->vertices(), nullptr);

	// Communities pivot vertices
	std::vector<std::vector<V>> commPivots(communities.size());

	// For each detected community extract related subgraph and vertex information 
	// for topological classes construction
	for (int i = 0; i < commVertexInfos.size(); i++)
	{
		std::shared_ptr<ISubGraph<V, W>> cluster = std::make_shared<SubGraph<V, W>>(communities[i]->all(), graph);

		_ce.evaluateCluster(commVertexInfos, cluster);

		commPivots[i] = _ps.selectPivots(commVertexInfos[i]);
	}

	std::valarray<W> intraClusterBC(globalBC);

	// Compute pivot contribution
	for (auto& community : commPivots)
	{
		for (auto& pivot : community)
		{
			std::valarray<W> pivotDependency = _ssb.singleSourceBrandes(pivot, graph);
			// TODO: Compute (d - communityDependency) * pivotClassCardinality
			// TODO: Sum contribution to globalBC
		}
	}

	return globalBC;
}

#endif