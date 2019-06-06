#ifndef FASTBC_KMEANS_IKMEANS
#define FASTBC_KMEANS_IKMEANS

#include <brandes/IVertexInfo.h>

#include <functional>
#include <unordered_map>
#include <vector>

namespace fastbc {
	namespace kmeans {

		template<typename V>
		class IKMeans
		{
		public:

			/**
			 *	@brief Compute k centroids from given vertex map
			 * 
			 *	@details Vertex info map should contain information useful for the vertex
			 *			 comparator function to correctly compute vertex distance
			 * 
			 *	@param k Number of centroids to compute
			 *	@param vertexInfo Map of vertex indices and related vertex information
			 *	@return std::vector<V> Vector of k centroids (vertex indices from vertexInfo input)
			 */
			virtual std::vector<V> computeCentroids(
				int k,
				const std::unordered_map<V, std::shared_ptr<IVertexInfo>> vertexInfo) = 0;
		};

	}
}

#endif