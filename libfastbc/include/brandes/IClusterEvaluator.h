#ifndef FASTBC_BRANDES_ICLUSTEREVALUATOR_H
#define FASTBC_BRANDES_ICLUSTEREVALUATOR_H

#include <ISubGraph.h>
#include "VertexInfo.h"

#include <memory>
#include <vector>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class IClusterEvaluator
		{
		public:

			/**
			 *	@brief Evaluate given sub-graph computing internal exact BC and other vertices information
			 * 
			 *	@details Betweenness centrality for each vertex is calculated along with 
			 *			 information about distance from border vertices and number of 
			 *			 shortest paths through them
			 * 
			 *	@note clusterBC and globalVI must be already initialized with correct 
			 *		  size of the global graph referenced by cluster sub-graph
			 * 
			 *	@param clusterBC Computed BC value will be summed to given reference
			 *	@param globalVI A new VertexInfo will be allocated for each of sub-graph vertices
			 *	@param cluster Sub-graph to apply computation to
			 */
			virtual void evaluateCluster(
				std::vector<W>& clusterBC,
				std::vector<std::shared_ptr<VertexInfo<V, W>>>& globalVI,
				std::shared_ptr<const ISubGraph<V,W>> cluster) = 0;
		};

	}
}

#endif