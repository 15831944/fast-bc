#ifndef FASTBC_BRANDES_ICLUSTEREVALUATOR_H
#define FASTBC_BRANDES_ICLUSTEREVALUATOR_H

#include <ISubGraph.h>
#include "IVertexInfo.h"

#include <map>
#include <memory>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class IClusterEvaluator
		{
		public:

			/**
			 *	@brief Evaluate given cluster computing exact betweenness centrality on each vertex
			 *
			 *	@details Cluster evaluation will compute exact BC along with some 
			 *			 additional information about shortest paths to each of the
			 *			 border vertices
			 * 
			 *	@param cluster Sub-graph object to evaluate
			 *	@return std::map<V, std::shared_ptr<IVertexInfo>> Computed information for each cluster vertex
			 */
			virtual std::map<V, std::shared_ptr<IVertexInfo>> evaluateCluster(
				std::shared_ptr<const ISubGraph<V,W>> cluster) = 0;
		};

	}
}

#endif