#ifndef FASTBC_LOUVAIN_ICOMMUNITY_H
#define FASTBC_LOUVAIN_ICOMMUNITY_H

#include <IGraph.h>

#include <memory>

namespace fastbc {

	namespace louvain {

		template<typename V, typename W>
		class ICommunity
		{
		public:

			/**
			 *	@brief Get all vertex contained in this community
			 *
			 *	@return std::vector<V> Full list of vertex in this community
			 */
			virtual std::vector<V> all() const = 0;

			/**
			 *	@brief Check if given vertex is contained in this community
			 * 
			 *	@param vertex Vertex's index
			 *	@return bool True if vertex is contained in this community, false else
			 */
			virtual bool contains(V vertex) const = 0;

			/**
			 *	@brief Add given vertex to this community
			 *
			 *	@param vertex Vertex's index
			 */
			virtual void add(V vertex) = 0;

			/**
			 *	@brief Remove given vertex from this community
			 *
			 *	@param vertex Vertex's index
			 */
			virtual void remove(V vertex) = 0;

			/**
			 *	@brief Get reference graph where this community lives
			 *
			 *	@return std::shared_ptr<const IGraph> Complete graph where this community is located
			 */
			virtual std::shared_ptr<const IGraph<V, W>> referenceGraph() const = 0;

			virtual int size() const = 0;
		};

	}

}

#endif
