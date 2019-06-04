#ifndef FASTBC_BRANDES_ISSBRANDESBC_H
#define FASTBC_BRANDES_ISSBRANDESBC_H

#include <IGraph.h>

#include <memory>
#include <valarray>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class ISSBrandesBC
		{
		public:

			/**
			 *	@brief Compute exact partial betweenness centrality values from given source vertex
			 *
			 *	@param source Source vertex
			 *	@param graph Full graph object
			 *	@return std::valarray<W> Partial betweenness centrality value for each graph vertex
			 */
			virtual std::valarray<W> singleSourceBrandes(
				V source, 
				std::shared_ptr<const IGraph<V, W>> graph) = 0;
		};

	}
}

#endif