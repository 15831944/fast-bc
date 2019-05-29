#ifndef FASTBC_BRANDES_ICLUSTEREVALUATOR_H
#define FASTBC_BRANDES_ICLUSTEREVALUATOR_H

#include <IGraph.h>
#include "IVertexInfo.h"

#include <map>
#include <vector>

namespace fastbc {
	namespace brandes {

		class IClusterEvaluator
		{
		public:

			/**
			 *	@brief Evaluate given cluster computing exact betweenness centrality on each vertex
			 *
			 *	@details Cluster evaluation will compute exact BC along with some 
			 *			additional information about shortest paths to each of the
			 *			border vertices
			 *
			 *	@param V Type of vertex index
			 *	@param W Type of edge weight value
			 *	@param graph Full graph descriptor
			 *	@param cluster Vertices to evaluate as a cluster (sub-graph)
			 *	@param borderVertices Border vertices of given cluster
			 *	@return std::map<V, IVertexInfo> Computed information for each cluster vertex
			 */
			template<typename V, typename W>
			virtual std::map<V, IVertexInfo> evaluateCluster(
				std::shared_ptr<const IGraph<V,W>> graph, 
				const std::vector<V>& cluster, 
				const std::vector<V>& borderVertices) = 0;
		};

	}
}

#endif