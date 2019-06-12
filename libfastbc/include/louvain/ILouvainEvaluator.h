#ifndef FASTBC_LOUVAIN_ILOUVAINEVALUATOR_H
#define FASTBC_LOUVAIN_ILOUVAINEVALUATOR_H

#include <IDegreeGraph.h>
#include <louvain/Partition.h>

namespace fastbc {
	namespace louvain {

		template<typename V, typename W>
		class ILouvainEvaluator
		{
		public:

			/**
			 *	@brief Evaluate given graph to create vertices communities using Louvain algorithm
			 *
			 *	@param graph Graph to evaluate
			 *	@return std::vector<std::vector<V>> Vertices communities computed from given graph
			 */				
			virtual std::vector<std::vector<V>> evaluateGraph(std::shared_ptr<const IDegreeGraph<V,W>> graph) = 0;
		};

	}
}

#endif
