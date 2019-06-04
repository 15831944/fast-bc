#ifndef FASTBC_ISUBGRAPH_H
#define FASTBC_ISUBGRAPH_H

#include <IGraph.h>

#include <memory>
#include <set>

namespace fastbc {

	template<typename V, typename W>
	class ISubGraph : public IGraph<V,W>
	{
	public:

		/**
		 *	@brief Get list of border vertices of this sub-graph
		 * 
		 *	@return Border vertices' indices iterators
		 */
		virtual const std::set<V>& borders() const = 0;

		/**
		 *	@brief Check if given vertex is at the border of this sub-graph
		 * 
		 *	@note A border vertex is a vertex such that it has at least one outgoing edge
		 *		  connected to a vertex outside the sub-graph.
		 * 
		 *	@param vertex Vertex index
		 * 
		 *	@return True if given vertex is at border, false else
		 */
		virtual bool isBorder(V vertex) const = 0;

		/**
		 *	@brief Get full graph where this sub-graph lies in
		 * 
		 *	@return Complete graph
		 */
		virtual std::shared_ptr<const IGraph<V, W>> referenceGraph() const = 0;
	};

}

#endif 