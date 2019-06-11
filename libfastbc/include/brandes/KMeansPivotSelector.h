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
			KMeansPivotSelector(
				std::shared_ptr<IPivotSelector<V, W>> exactPivotSelector,
				std::shared_ptr<kmeans::IKMeans<V, W>> kmeans,
				double kFrac,
				W stopVariance = 0,
				size_t maxIteration = 100);

			std::pair<std::vector<V>, std::vector<V>> selectPivots(
				const std::valarray<W>& globalBC,
				const std::vector<std::shared_ptr<VertexInfo<V, W>>>& verticesInfo,
				const std::set<V>& vertices,
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
	const std::valarray<W>& globalBC,
	const std::vector<std::shared_ptr<VertexInfo<V, W>>>& verticesInfo,
	const std::set<V>& vertices,
	const std::set<V>& borders)
{
	const auto[pivotIndexCluster, pivotClassCluster] = 
		_exactPS->selectPivots(globalBC, verticesInfo, vertices, borders);

	int k = std::max((int)(pivotIndexCluster.size() * _kFrac), 1);

	SPDLOG_INFO("Aggregating {} pivots in {} super-classes", 
		pivotIndexCluster.size(), k);

	std::pair<std::vector<V>, std::vector<V>> pivotWeight = 
		_kmeans->computeCentroids(k, pivotIndexCluster, pivotClassCluster, verticesInfo, 
			_stopVariance, _maxIteration);

#ifndef FASTBC_BRANDES_KMENS_PIVOT_ALLOW_DUPLICATED
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
		SPDLOG_ERROR("Removed {} duplicated pivots from current cluster", 
			duplicates);
	}
#endif
	

	return pivotWeight;
}

#endif