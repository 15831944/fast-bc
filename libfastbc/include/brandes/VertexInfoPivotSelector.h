#ifndef FASTBC_BRANDES_VERTEXINFOPIVOTSELECTOR_H
#define FASTBC_BRANDES_VERTEXINFOPIVOTSELECTOR_H

#include "IPivotSelector.h"

#include <memory>
#include <spdlog/spdlog.h>
#include <vector>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class VertexInfoPivotSelector : public IPivotSelector<V, W>
		{
		public:

			std::pair<std::vector<V>, std::vector<V>> selectPivots(
				const std::valarray<W>& globalBC,
				const std::vector<std::shared_ptr<VertexInfo<V, W>>>& verticesInfo,
				const std::set<V>& vertices,
				const std::set<V>& borders) override;
		};

	}
}

template<typename V, typename W>
std::pair<std::vector<V>, std::vector<V>> fastbc::brandes::VertexInfoPivotSelector<V, W>::selectPivots(
	const std::valarray<W>& globalBC,
	const std::vector<std::shared_ptr<VertexInfo<V, W>>>& verticesInfo,
	const std::set<V>& vertices,
	const std::set<V>& borders)
{
	// Vertex info class representative
	std::vector<std::shared_ptr<VertexInfo<V, W>>> classes;

	// Vertices for each class
	std::vector<std::vector<V>> classMembers;

	// Add each vertex to a class
	for (const auto& v : vertices)
	{
		std::shared_ptr<VertexInfo<V, W>> vVI = verticesInfo[v];
		
		// Normalize VI before comparison to allow correct class aggregation
		vVI->normalize();

		// Check if a suitable class already exists
		bool classExists = false;
		for (V ci = 0; ci < classes.size(); ++ci)
		{
			const auto& cVI = classes[ci];

			if (*cVI == *vVI)
			{
				classMembers[ci].push_back(v);
				classExists = true;
				break;
			}
		}

		// If no class exists for current vertex generate a new one
		if (!classExists)
		{
			classes.push_back(vVI);
			classMembers.push_back(std::vector<V>({ v }));
		}
	}

	SPDLOG_INFO("Found {} topological classes in current cluster", classes.size());

	// Classes pivot and cardinality
	std::pair<std::vector<V>, std::vector<V>> pivot;

	// Update each vertex class cardinality and select vertex with minimum BC as pivot
	for (const auto& classM : classMembers)
	{
		// Try find first non-border node
		int j = 0;
		while (j < classM.size())
		{
			if (borders.find(classM[j]) != borders.end())
			{
				++j;
			}
			else
			{
				break;
			}
		}

		// If class has only border vertices
		if (j == classM.size())
		{
			
#ifdef FASTBC_BRANDES_ENABLE_PIVOT_BORDER
			SPDLOG_WARN("Topological class contains only border vertices: selecting first as pivot");
			pivot.first.push_back(classM[0]);
			pivot.second.push_back(classM.size());
#else
			SPDLOG_WARN("Topological class contains only border vertices: no pivot was selected");
#endif
			continue;
		}

		// First eligible pivot node
		V minV = classM[j];

		for (; j < classM.size(); ++j)
		{
			V v = classM[j];

			// ONLY NON-BORDER NODES CAN BE SELECTED AS PIVOTS
			if (borders.find(v) == borders.end() && globalBC[v] < globalBC[minV])
			{
				minV = v;
			}
		}

		// Store selected pivot and related class cardinality
		pivot.first.push_back(minV);
		pivot.second.push_back(classM.size());
	}

	return pivot;
}

#endif