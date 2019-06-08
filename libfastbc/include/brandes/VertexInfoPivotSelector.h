#ifndef FASTBC_BRANDES_VERTEXINFOPIVOTSELECTOR_H
#define FASTBC_BRANDES_VERTEXINFOPIVOTSELECTOR_H

#include "IPivotSelector.h"

#include <memory>
#include <vector>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class VertexInfoPivotSelector : public IPivotSelector<V, W>
		{
		public:

			std::vector<V> selectPivots(
				const std::valarray<W>& globalBC,
				const std::vector<std::shared_ptr<VertexInfo<V, W>>>& verticesInfo,
				std::valarray<W>& verticesClassCardinality,
				const std::set<V>& vertices,
				const std::set<V>& borders) override;
		};

	}
}

template<typename V, typename W>
std::vector<V> fastbc::brandes::VertexInfoPivotSelector<V, W>::selectPivots(
	const std::valarray<W>& globalBC,
	const std::vector<std::shared_ptr<VertexInfo<V, W>>>& verticesInfo,
	std::valarray<W>& verticesClassCardinality,
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

	// Classes pivot result
	std::vector<V> pivotVertices(classes.size());

	// Update each vertex class cardinality and select vertex with minimum BC as pivot
	for (V i = 0; i < classes.size(); ++i)
	{
		const auto& classM = classMembers[i];

		// Find first non border node
		int j = 0;
		while(borders.find(classM[j]) != borders.end()) j++;
		V minV;
		if(j < classM.size())
			minV = classM[j];
		else
			minV = classM[0];

		for (auto& v : classM)
		{
			verticesClassCardinality[v] = classM.size();

			// ONLY NON-BORDER NODES CAN BE SELECTED AS PIVOTS
			if (borders.find(v) == borders.end() && globalBC[v] < globalBC[minV])
			{
				minV = v;
			}
		}

		pivotVertices[i] = minV;
	}

	return pivotVertices;
}

#endif