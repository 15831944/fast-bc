#ifndef FASTBC_BRANDES_IVERTEXINFO_H
#define FASTBC_BRANDES_IVERTEXINFO_H

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class IVertexInfo
		{
		public:
			/**
			 *	@brief Get current centrality value
			 *
			 *	@return W Vertex centrality value
			 */
			virtual W getCentrality() const = 0;

			/**
			 *	@brief Set new centrality value
			 *
			 *	@param c New centrality value 
			 */
			virtual void setCentrality(W c) = 0;

			/**
			 *	@brief Get length for shortest path to border vertex
			 *
			 *	@param borderVertex Border vertex destination of required shortest path
			 *	@return W Shortest path length
			 */
			virtual W distanceFromBorder(V borderVertex) const = 0;

			/**
			 *	@brief Get shortest path count to given border vertex
			 *
			 *	@param borderVertex Border vertex destination of required shortest path count
			 *	@return V Number of shortest path going to border vertex
			 */
			virtual V spCountToBorder(V borderVertex) const = 0;

			/**
			 *	@brief Compute topological distance between this and another vertex information
			 *	
			 *	@details Topological distance is represented by the Euclidean distance between 
			 *			 the vector of border distances/shortest path count of this and vi 
			 *			 vertices information
			 * 
			 *	@param vi Vertex info to compare to this
			 *	@return int Zero if topological charateristics are equal, distance else
			 */
			virtual int topologicalDistance(const IVertexInfo<V,W>& vi) const = 0;

			static int topologicalDistance(const IVertexInfo<V,W>& a, const IVertexInfo& b);
		};

	}
}

template<typename V, typename W>
int fastbc::brandes::IVertextInfo<V, W>::topologicalDistance(const IVertexInfo<V,W>& a, const IVertexInfo<V,W>& b)
{
	return a.topologicalDistance(b);
}

#endif // !FASTBC_BRANDES_IVERTEXINFO_H
