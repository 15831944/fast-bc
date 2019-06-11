#ifndef FASTBC_KMEANS_IKMEANS
#define FASTBC_KMEANS_IKMEANS

#include <brandes/VertexInfo.h>

#include <memory>
#include <vector>
#include <utility>

namespace fastbc {
	namespace kmeans {

		template<typename V, typename W>
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
			 *	@param vertices Vertex indices to compute clusters on
			 *	@param weights Vertex weights to consider during clusters computation
			 *	@param vertexInfo Map of vertex indices and related vertex information
			 *	@param minVariance Minimum subsequent iteration centroids variance to consider
			 *	@param maxIteration Maximum number of iterations allowed
			 *	@return std::pair<std::vector<V>, std::vector<V>> Vector of k centroids and related weights
			 */
			virtual std::pair<std::vector<V>, std::vector<V>> computeCentroids(
				int k,
				const std::vector<V>& vertices,
				const std::vector<V>& weights,
				const std::vector<std::shared_ptr<brandes::VertexInfo<V, W>>>& vertexInfo,
				W minVariance = 0,
				size_t maxIteration = 100) = 0;
		};

	}
}

#endif