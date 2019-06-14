#ifndef FASTBC_BRANDES_IBRANDESBC_H
#define FASTBC_BRANDES_IBRANDESBC_H

#include <IGraph.h>

#include <memory>
#include <vector>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class IBrandesBC
		{
		public:

			/**
			 * 	@brief Compute betweenness centrality of graph using Brandes' algorithm as main routine
			 * 
			 * 	@note graph must be a compelte graph: vertex indices from 0 to graph->vertices().size()
			 * 
			 * 	@param graph Complete graph to compute BC for
			 */
			virtual std::vector<W> computeBC(const std::shared_ptr<const IGraph<V, W>> graph) = 0;
		};

	}
}

#endif