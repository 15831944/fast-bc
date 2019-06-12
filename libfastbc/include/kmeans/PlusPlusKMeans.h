#ifndef FASTBC_KMEANS_PLUSPLUSKMEANS_H
#define FASTBC_KMEANS_PLUSPLUSKMEANS_H

#include "IKMeans.h"

namespace fastbc {
	namespace kmeans {

		template<typename V, typename W>
		class PlusPlusKMeans : public IKMeans<V, W>
		{
		public:
			std::pair<std::vector<V>, std::vector<V>> computeCentroids(
				int k,
				const std::vector<V>& vertices,
				const std::vector<V>& weights,
				const std::vector<std::shared_ptr<brandes::VertexInfo<V, W>>>& vertexInfo,
				W minVariance = 0,
				size_t maxIteration = 100) override;

		private:

			std::vector<V> _initPlusPlus(
				int k,
				const std::vector<V>& vertices,
				const std::vector<std::shared_ptr<brandes::VertexInfo<V, W>>>& vertexInfo);

			W _centroidVariance(
				const std::vector<V>& oldCentroid,
				const std::vector<V>& newCentroid,
				const std::vector<std::shared_ptr<brandes::VertexInfo<V, W>>>& vertexInfo);

			struct InfoCluster { 
				brandes::VertexInfo<W, W> centroidInfo;
				std::vector<V> vIndices;

				InfoCluster(int borderCount) : centroidInfo(borderCount) {}

				InfoCluster& operator+=(const InfoCluster& other)
				{
					centroidInfo += other.centroidInfo;
					vIndices.insert(vIndices.end(), other.vIndices.begin(), other.vIndices.end());
					return *this;
				}

				void reset()
				{
					centroidInfo.reset();
					vIndices.clear();
				}
			};
			#pragma omp declare reduction(+: InfoCluster: omp_out += omp_in) \
				initializer(omp_priv = InfoCluster(omp_orig.centroidInfo.borders()))

			struct VertexDistance { 
				V vertex; 
				W distance; 

				VertexDistance(V v, W d) : vertex(v), distance(d) {}
			};
			#pragma omp declare reduction(min: VertexDistance: \
				omp_out = omp_out.distance < omp_in.distance ? omp_out : omp_in) \
				initializer(omp_priv = VertexDistance(omp_orig.vertex, omp_orig.distance))

		};

	}
}

template<typename V, typename W>
std::pair<std::vector<V>, std::vector<V>>
fastbc::kmeans::PlusPlusKMeans<V, W>::computeCentroids(
	int k,
	const std::vector<V>& vertices,
	const std::vector<V>& weights,
	const std::vector<std::shared_ptr<brandes::VertexInfo<V, W>>>& vertexInfo,
	W minVariance,
	size_t maxIteration)
{
	// Current centroids vector
	std::vector<V> newCentroid = _initPlusPlus(k, vertices, vertexInfo);
	std::vector<V> centroid(newCentroid.size());

	std::vector<InfoCluster> infoCluster(centroid.size(), InfoCluster(vertexInfo[vertices[0]]->borders()));
	InfoCluster* _infoCluster = infoCluster.data();
	size_t _infoClusterSize = infoCluster.size();

	size_t iteration = 0;
	do {
		++iteration;
		centroid = newCentroid;
		for (auto& ic : infoCluster) { ic.reset(); }

		// Associate each vertex to nearest cluster
		#pragma omp parallel for reduction(+:_infoCluster[:_infoClusterSize])
		for (size_t v = 0; v < vertices.size(); ++v)
		{
			struct VertexDistance minC(0, 
				vertexInfo[centroid[0]]->contributionDistance(*vertexInfo[vertices[v]]));

			// Select nearest cluster to current vertex
			#pragma omp simd reduction(min:minC)
			for (int c = 1; c < centroid.size(); ++c)
			{
				W dist = vertexInfo[centroid[c]]->contributionDistance(*vertexInfo[vertices[v]]);

				if (dist < minC.distance)
				{
					minC.vertex = c;
					minC.distance = dist;
				}
			}

			// Store vertex association to selected cluster
			_infoCluster[minC.vertex].centroidInfo += *vertexInfo[vertices[v]];
			_infoCluster[minC.vertex].vIndices.push_back(v);
		}

		// Choose new centroids for each computed cluster
		for (int c = 0; c < centroid.size(); ++c)
		{
			auto& ic = infoCluster[c];
			if (ic.vIndices.empty())
			{
				continue;
			}

			ic.centroidInfo /= ic.vIndices.size();

			struct VertexDistance minV(ic.vIndices[0],
				ic.centroidInfo.contributionDistance(*vertexInfo[vertices[ic.vIndices[0]]]));

			// New centroid will be the nearest existing vertex to computed centroid
			#pragma omp simd reduction(min:minV)
			for (size_t v = 1; v < ic.vIndices.size(); ++v)
			{
				W dist = ic.centroidInfo.contributionDistance(*vertexInfo[vertices[ic.vIndices[v]]]);

				if (dist < minV.distance)
				{
					minV.vertex = ic.vIndices[v];
					minV.distance = dist;
				}
			}

			newCentroid[c] = vertices[minV.vertex];
		}

	// Iterate until no significant change in centroids is detected or max iteration is reached
	} while (_centroidVariance(centroid, newCentroid, vertexInfo) > minVariance
		&& iteration <= maxIteration);


	std::pair<std::vector<V>, std::vector<V>> centroidWeights = 
		std::make_pair(centroid, std::vector<V>(centroid.size(), 0));

	// Compute each centroid final weight
	V* cWeights = centroidWeights.second.data();
	size_t cWeightsSize = centroidWeights.second.size();
	for (int c = 0; c < centroid.size(); ++c)
	{
		#pragma omp simd reduction(+:cWeights[:cWeightsSize])
		for (size_t v = 0; v < infoCluster[c].vIndices.size(); ++v)
		{
			cWeights[c] += weights[infoCluster[c].vIndices[v]];
		}
	}

	return centroidWeights;
}

template<typename V, typename W>
std::vector<V>
fastbc::kmeans::PlusPlusKMeans<V, W>::_initPlusPlus(
	int k,
	const std::vector<V>& vertices,
	const std::vector<std::shared_ptr<brandes::VertexInfo<V, W>>>& vertexInfo)
{
	std::vector<V> centroid(k);
	centroid[0] = vertices[0];

	std::vector<W> cDist(vertices.size(), 0);
	for (int i = 1; i < k; ++i)
	{
		std::shared_ptr<brandes::VertexInfo<V, W>> lastCentroid = vertexInfo[centroid[i - 1]];
		double p = 1.0 / i;
		double _p = 1.0 - p;

		V farthestV = 0;

		for (int v = 0; v < vertices.size(); ++v)
		{
			// Update distance from prevoiusly selected centroids
			cDist[v] = cDist[v] * _p + lastCentroid->contributionDistance(*vertexInfo[vertices[v]]) * p;

			// Update farthest from existing centroids
			if (cDist[v] > cDist[farthestV])
			{
				farthestV = v;
			}
		}

		centroid[i] = vertices[farthestV];
	}

	return centroid;
}

template<typename V, typename W>
W fastbc::kmeans::PlusPlusKMeans<V, W>::_centroidVariance(
	const std::vector<V>& oldCentroid,
	const std::vector<V>& newCentroid,
	const std::vector<std::shared_ptr<brandes::VertexInfo<V, W>>>& vertexInfo)
{
	W maxVariance = 0;

	#pragma omp simd reduction(max:maxVariance)
	for (int c = 0; c < oldCentroid.size(); ++c)
	{
		W variance = vertexInfo[oldCentroid[c]]->contributionDistance(*vertexInfo[newCentroid[c]]);

		if (variance > maxVariance)
		{
			maxVariance = variance;
		}
	}

	return maxVariance;
}

#endif