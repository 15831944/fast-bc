#ifndef FASTBC_BRANDES_ISSBRANDESBC_H
#define FASTBC_BRANDES_ISSBRANDESBC_H

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class ISSBrandesBC
		{
		public:

			/**
			 *	@brief Compute exact partial betweenness centrality values from given source vertex
			 *
			 *	@param vertex Source vertex
			 *	@param graph Full graph object
			 *	@return std::vector<W> Partial betweenness centrality value for each graph vertex
			 */
			virtual std::vector<W> singleSourceBrandes(
				V source, 
				std::shared_ptr<const IGraph<V, W>> graph) = 0;
		};

	}
}

#endif