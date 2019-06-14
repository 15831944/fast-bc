#ifndef FASTBC_IGRAPHPARTITION_H
#define FASTBC_IGRAPHPARTITION_H

#include "IDegreeGraph.h"

namespace fastbc {

	template<typename V, typename W>
	class IGraphPartition
	{
	public:
	
		/**
		 *	@brief Evaluate given graph to create vertices clusters
		 *
		 *	@param graph Graph to evaluate
		 *	@return std::vector<std::vector<V>> Vertices communities computed from given graph
		 */				
		virtual std::vector<std::vector<V>> partitionGraph(std::shared_ptr<const IDegreeGraph<V,W>> graph) = 0;
	};
}

#endif
