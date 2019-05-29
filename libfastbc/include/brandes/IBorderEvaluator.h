#ifndef FASTBC_BRANDES_IBORDEREVALUATOR_H
#define FASTBC_BRANDES_IBORDEREVALUATOR_H


namespace fastbc {
	namespace brandes {

		class IBorderEvaluator
		{
		public:

			/**
			 *	@brief Find border vertices in given sub-graph
			 *
			 *	@param V Type of vertex index
			 *	@param W Type of edge weight value
			 *	@param graph Graph instance
			 *	@param subGraph Vertex indices to consider inside of sub-graph
			 *	@return std::vector<V> Indices of border vertices
			 */
			template<typename V, typename W>
			virtual std::vector<V> findBorders(
				std::shared_ptr<const IGraph<V, W>> graph,
				const std::vector<V>& subGraph) = 0;
		};

	}
}

#endif