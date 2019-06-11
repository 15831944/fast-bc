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

	std::vector<std::vector<V>> clusterVerticesIndex(centroid.size());

	size_t iteration = 0;

	do {
		++iteration;
		centroid = newCentroid;
		for (auto& cvi : clusterVerticesIndex) { cvi.clear(); }
		
		std::vector<brandes::VertexInfo<W, W>> centroidInfo(centroid.size(), brandes::VertexInfo<W, W>(vertexInfo[vertices[0]]->borders()));

		// Associate each vertex to nearest cluster
		for (size_t v = 0; v < vertices.size(); ++v)
		{
			int nearestC = 0;
			W nearestDist = vertexInfo[centroid[nearestC]]->contributionDistance(*vertexInfo[vertices[v]]);

			// Select nearest cluster to current vertex
			for (int c = 1; c < centroid.size(); ++c)
			{
				W dist = vertexInfo[centroid[c]]->contributionDistance(*vertexInfo[vertices[v]]);

				if (dist < nearestDist)
				{
					nearestC = c;
					nearestDist = dist;
				}
			}

			// Store vertex association to selected cluster
			clusterVerticesIndex[nearestC].push_back(v);
			centroidInfo[nearestC] += *vertexInfo[vertices[v]];
		}

		// Choose new centroids for each computed cluster
		for (int c = 0; c < centroid.size(); ++c)
		{
			if (clusterVerticesIndex[c].empty())
			{
				continue;
			}

			centroidInfo[c] /= clusterVerticesIndex[c].size();

			V nearestV = clusterVerticesIndex[c][0];
			W nearestDist = centroidInfo[c].contributionDistance(*vertexInfo[vertices[nearestV]]);

			// New centroid will be the nearest existing vertex to computed centroid
			for (size_t v = 1; v < clusterVerticesIndex[c].size(); ++v)
			{
				W dist = centroidInfo[c].contributionDistance(*vertexInfo[vertices[clusterVerticesIndex[c][v]]]);

				if (dist < nearestDist)
				{
					nearestV = clusterVerticesIndex[c][v];
					nearestDist = dist;
				}
			}

			newCentroid[c] = vertices[nearestV];
		}

	// Iterate until no significant change in centroids is detected or max iteration is reached
	} while (_centroidVariance(centroid, newCentroid, vertexInfo) > minVariance
		&& iteration <= maxIteration);


	std::pair<std::vector<V>, std::vector<V>> centroidWeights = 
		std::make_pair(centroid, std::vector<V>(centroid.size(), 0));

	// Compute each centroid final weight
	for (int c = 0; c < centroid.size(); ++c)
	{
		for (size_t v = 0; v < clusterVerticesIndex[c].size(); ++v)
		{
			centroidWeights.second[c] += weights[clusterVerticesIndex[c][v]];
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