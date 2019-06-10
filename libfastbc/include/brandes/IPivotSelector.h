#ifndef FASTBC_BRANDES_IPIVOTSELECTOR_H
#define FASTBC_BRANDES_IPIVOTSELECTOR_H

#include "VertexInfo.h"

#include <memory>
#include <set>
#include <valarray>
#include <vector>
#include <utility>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class IPivotSelector
		{
		public:

			/**
			 *	@brief Select pivot for topological classes based on given vertices information
			 * 
			 *	@details Generated pivots are vertices with smallest BC in their class and not
			 *			 border; a class is composed of vertices with equal vertex information 
			 * 
			 *	@note Given VertexInfo references will be normalized and class caradinality will
			 *		  be updated with correct value during the call
			 * 
			 *	@param globalBC Betweenness centrality value for each vertex
			 *	@param verticesInfo Vertex information for each vertex
			 *	@param vertices Vertices to be considered in the computation
			 *	@param borders Vertices not to be considered as pivot
			 *	@return std::pair<std::vector<V>, std::vector<V>> Selected pivot vertex indices and related class cardinality
			 */
			virtual std::pair<std::vector<V>, std::vector<V>> selectPivots(
				const std::valarray<W>& globalBC, 
				const std::vector<std::shared_ptr<VertexInfo<V, W>>>& verticesInfo,
				const std::set<V>& vertices,
				const std::set<V>& borders) = 0;
		};

	}
}

#endif // !FASTBC_BRANDES_IPIVOTSELECTOR_H
