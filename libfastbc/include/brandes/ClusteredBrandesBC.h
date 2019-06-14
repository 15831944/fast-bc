#ifndef FASTBC_BRANDES_CLUSTEREDBRANDESBC_H
#define FASTBC_BRANDES_CLUSTEREDBRANDESBC_H

#include "IBrandesBC.h"
#include "IClusterEvaluator.h"
#include "ISSBrandesBC.h"
#include "IPivotSelector.h"
#include "VertexInfo.h"
#include <IGraphPartition.h>
#include <SubGraph.h>

#include <memory>
#include <spdlog/spdlog.h>
#include <vector>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class ClusteredBrandeBC : public IBrandesBC<V, W>
		{
		public:
			/*
			 *	@brief Initialize a clustered brandes BC computer
			 * 
			 * 	@details The object will perform a clustered Brandes' BC computation
			 * 			 using given algorithm instances for clusters creation and
			 * 			 evaluation and BC computation
			 * 
			 * 	@param gp Graph partition creator
			 * 	@param ce Cluster BC evaluator
			 * 	@param ssb Single source Brandes' BC computer
			 * 	@param ps Pivot selector to use on computed clusters
			 */
			ClusteredBrandeBC(
				std::shared_ptr<IGraphPartition<V, W>> gp,
				std::shared_ptr<IClusterEvaluator<V, W>> ce,
				std::shared_ptr<ISSBrandesBC<V, W>> ssb,
				std::shared_ptr<IPivotSelector<V, W>> ps);

			std::vector<W> computeBC(const std::shared_ptr<const IGraph<V, W>> graph) override;

		private:
			std::shared_ptr<IGraphPartition<V, W>> _gp;
			std::shared_ptr<IClusterEvaluator<V, W>> _ce;
			std::shared_ptr<ISSBrandesBC<V, W>> _ssb;
			std::shared_ptr<IPivotSelector<V, W>> _ps;
		};

	}
}

template<typename V, typename W>
fastbc::brandes::ClusteredBrandeBC<V, W>::ClusteredBrandeBC(
	std::shared_ptr<fastbc::IGraphPartition<V, W>> gp,
	std::shared_ptr<fastbc::brandes::IClusterEvaluator<V, W>> ce,
	std::shared_ptr<fastbc::brandes::ISSBrandesBC<V, W>> ssb,
	std::shared_ptr<fastbc::brandes::IPivotSelector<V, W>> ps)
	: _gp(gp), _ce(ce), _ssb(ssb), _ps(ps)
{
}

template<typename V, typename W>
std::vector<W> fastbc::brandes::ClusteredBrandeBC<V, W>::computeBC(
	const std::shared_ptr<const fastbc::IGraph<V, W>> graph)
{
	// Global betweenness centrality storage
	std::vector<W> globalBC(graph->vertices().size(), (W)0);

	// Vertices topological information about their own cluster border vertices
	std::vector<std::shared_ptr<VertexInfo<V, W>>> verticesInfo(graph->vertices().size(), nullptr);

	// Computed subgraph and border vertices from each vertices community
	std::vector<std::shared_ptr<ISubGraph<V, W>>> cluster;

	// Pivot vertices and related class cardinality for each cluster
	std::vector<std::pair<std::vector<V>, std::vector<V>>> pivotsCluster;

	// Compute graph partition using Louvain communities detection algorithm
	SPDLOG_INFO("Computing clusters with Louvain algorithm...");
	std::vector<std::vector<V>> communities = 
		_gp->partitionGraph(std::static_pointer_cast<const IDegreeGraph<V, W>>(graph));

	SPDLOG_INFO("Graph partitioned in {} clusters", communities.size());
	cluster.resize(communities.size());
	pivotsCluster.resize(communities.size());

	// For each detected community compute related sub-graph, evaluate it for internal BC
	// and perform topological analysis to get pivots and vertices class cardinality
	SPDLOG_INFO("Evaluating intra cluster BC...");
	#pragma omp parallel for
	for (int i = 0; i < cluster.size(); i++)
	{
		cluster[i] = std::make_shared<SubGraph<V, W>>(communities[i], graph);

		SPDLOG_DEBUG("Evaluating BC on cluster {}: {} vertices ({} borders), {} edges", 
			i, cluster[i]->vertices().size(), cluster[i]->borders().size(), cluster[i]->edges());
		
#ifndef FASTBC_BRANDES_CLUSTERED_IGNORE_UNCONNECTED
		if (cluster[i]->borders().empty())
		{
			SPDLOG_WARN("Cluster {} ({} vertices, {} edges) is disconnected from the rest of the graph.", 
				i, cluster[i]->vertices().size(), cluster[i]->edges());
		}
#else
		if (!cluster[i]->borders().empty())
		{
#endif
		
		_ce->evaluateCluster(globalBC, verticesInfo, cluster[i]);

		pivotsCluster[i] = _ps->selectPivots(
			globalBC, verticesInfo, 
			cluster[i]->vertices(), cluster[i]->borders());

		SPDLOG_DEBUG("Selected {} vertices as pivots in cluster {}", pivotsCluster[i].first.size(), i);
		
#ifdef FASTBC_BRANDES_CLUSTERED_IGNORE_UNCONNECTED
		}
#endif
	}

	// Store computed intra-cluster BC for corrections on 
	// following global BC computation step
	std::vector<W> intraClusterBC(globalBC);
	
	// Print total number of selected pivots
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_INFO
	size_t pivotCount = 0;
	#pragma omp simd reduction(+:pivotCount)
	for (size_t c = 0; c < pivotsCluster.size(); ++c)
	{
		pivotCount += pivotsCluster[c].first.size();
	}
	SPDLOG_INFO("Computing global BC from {} pivots...", pivotCount);
#endif

	// Compute global dependecy contribution for each selected pivot
	W* _globalBC = globalBC.data();
	size_t _globalBCsize = globalBC.size();
	for (size_t c = 0; c < cluster.size(); ++c)
	{
		#pragma omp parallel for reduction(+:_globalBC[:_globalBCsize])
		for (size_t p = 0; p < pivotsCluster[c].first.size(); ++p)
		{
			SPDLOG_DEBUG("Computing SSSP from pivot vertex {}", pivotsCluster[c].first[p]);
			std::vector<W> pivotDependency = 
				_ssb->singleSourceBrandes(pivotsCluster[c].first[p], graph);

			// Sum pivot dependecy to all vertices
			#pragma omp simd
			for(size_t v = 0; v < _globalBCsize; ++v)
			{
				_globalBC[v] += pivotDependency[v] * (W)(pivotsCluster[c].second[p]);
			}

			// Subtract duplicate dependency from current pivot's cluster vertices
			#pragma omp simd
			for (size_t vIndex = 0; vIndex < cluster[c]->vertices().size(); ++vIndex)
			{
				const V& v = cluster[c]->vertices()[vIndex];

				_globalBC[v] -= intraClusterBC[v] * (W)(pivotsCluster[c].second[p]);
			}
		}
	}

	return globalBC;
}

#endif