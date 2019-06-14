#ifndef FASTBC_BRANDES_KMEANSPIVOTSELECTOR_H
#define FASTBC_BRANDES_KMEANSPIVOTSELECTOR_H

#include "IPivotSelector.h"
#include <kmeans/IKMeans.h>

#include <algorithm>
#include <memory>
#include <spdlog/spdlog.h>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class KMeansPivotSelector : public IPivotSelector<V, W>
		{
		public:

			/**
			 * 	@brief Initialize a kmeans pivot selector
			 * 
			 * 	@details Pivot selection is first computed considering exact VertexInfo
			 * 			 topological classes, then kmeans is applied to select a subset
			 * 			 of those pivots based on VertexInfo::squaredDistance
			 * 
			 * 	@param exactPivotSelector Exact pivot selector used in first step
			 * 	@param kmeans KMeans computer used in second step
			 * 	@param kFrac Fraction of exact pivots to exctract using kmeans in second step
			 * 	@param stopVariance Minimum allowed variance between pivots set to trigger a new kmeans iteration
			 * 	@param maxIteration Maximumallowed kmeans iterations
			 */
			KMeansPivotSelector(
				std::shared_ptr<IPivotSelector<V, W>> exactPivotSelector,
				std::shared_ptr<kmeans::IKMeans<V, W>> kmeans,
				double kFrac,
				W stopVariance = 0,
				size_t maxIteration = 100);

			std::pair<std::vector<V>, std::vector<V>> selectPivots(
				const std::vector<W>& globalBC,
				const std::vector<std::shared_ptr<VertexInfo<V, W>>>& verticesInfo,
				const std::vector<V>& vertices,
				const std::set<V>& borders) override;

		private:
			std::shared_ptr<IPivotSelector<V, W>> _exactPS;
			std::shared_ptr<kmeans::IKMeans<V, W>> _kmeans;
			const double _kFrac;
			const W _stopVariance;
			const size_t _maxIteration;
		};
	}
}

template<typename V, typename W>
fastbc::brandes::KMeansPivotSelector<V, W>::KMeansPivotSelector(
	std::shared_ptr<IPivotSelector<V, W>> exactPivotSelector,
	std::shared_ptr<kmeans::IKMeans<V, W>> kmeans,
	double kFrac,
	W stopVariance,
	size_t maxIteration)
	: _exactPS(exactPivotSelector), 
	_kmeans(kmeans), 
	_kFrac(kFrac), 
	_stopVariance(stopVariance),
	_maxIteration(maxIteration)
{
	if (_kFrac < 0.0 || _kFrac > 1.0)
	{
		throw std::invalid_argument("Given KFrac parameter is out of bounds");
	}

	if(maxIteration < 100)
	{
		SPDLOG_WARN("Given max iteration for kmeans pivot selection is low ({})", _maxIteration);
	}
}

template<typename V, typename W>
std::pair<std::vector<V>, std::vector<V>> 
fastbc::brandes::KMeansPivotSelector<V, W>::selectPivots(
	const std::vector<W>& globalBC,
	const std::vector<std::shared_ptr<VertexInfo<V, W>>>& verticesInfo,
	const std::vector<V>& vertices,
	const std::set<V>& borders)
{
	// Compute exact topological classes and their pivots
	const auto[pivotIndexCluster, pivotClassCluster] = 
		_exactPS->selectPivots(globalBC, verticesInfo, vertices, borders);

	// Compute pivots subset cardinality
	int k = std::max((int)(pivotIndexCluster.size() * _kFrac), 1);

	SPDLOG_TRACE("Aggregating {} pivots in {} super-classes", 
		pivotIndexCluster.size(), k);

	// Compute pivots subset through kmeans algorithm
	// BE AWARE: duplicated pivots can result from kmeans due to the algorithm euristic nature
	std::pair<std::vector<V>, std::vector<V>> pivotWeight = 
		_kmeans->computeCentroids(k, pivotIndexCluster, pivotClassCluster, verticesInfo, 
			_stopVariance, _maxIteration);

#ifndef FASTBC_BRANDES_KMENS_PIVOT_ALLOW_DUPLICATED
	// Remove duplicated pivots from kmeans result
	std::set<V> uniquePivots;
	size_t duplicates = 0;
	auto pIT = pivotWeight.first.begin();
	auto wIT = pivotWeight.second.begin();
	while(pIT != pivotWeight.first.end())
	{
		if(uniquePivots.insert(*pIT).second)
		{
			++pIT;
			++wIT;
		}
		else
		{
			pIT = pivotWeight.first.erase(pIT);
			wIT = pivotWeight.second.erase(wIT);
			++duplicates;
		}
	}
	
	if(duplicates)
	{
		SPDLOG_TRACE("Removed {} duplicated pivots from current cluster", 
			duplicates);
	}
#endif
	
	return pivotWeight;
}

#endif