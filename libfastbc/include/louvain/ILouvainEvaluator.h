#ifndef FASTBC_LOUVAIN_ILOUVAINEVALUATOR_H
#define FASTBC_LOUVAIN_ILOUVAINEVALUATOR_H

#include <IGraph.h>
#include "ICommunity.h"

namespace fastbc {
	namespace louvain {

		class ILouvainEvaluator
		{
		public:

			/**
			 *	@brief Evaluate given graph to create vertices communities using Louvain algorithm
			 *
			 *	@param graph Graph to evaluate
			 *	@return std::vector<ICommunity<V,W>> Vertices communities computed from given graph
			 */
			template<typename V, typename W>
			virtual std::vector<ICommunity<V,W>> evaluateGraph(std::shared_ptr<const IGraph<V,W>> graph) = 0;
		};

	}
}

#endif