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
			 *	@brief Check if two vertex has same topological charactericìstics with respect to border nodes
			 *	
			 * @param vi Vertex info to compare to this
			 * @return bool True if vi has same topological characteristics, false else
			 */
			virtual bool compareTopologyClass(const IVertexInfo& vi) const = 0;
		};

	}
}

#endif // !FASTBC_BRANDES_IVERTEXINFO_H
