#include <catch2/catch.hpp>

#include <brandes/VertexInfoPivotSelector.h>

#include <brandes/VertexInfo.h>
#include <algorithm>
#include <memory>
#include <set>
#include <valarray>
#include <vector>

using namespace fastbc::brandes;

TEST_CASE("Pivot selection", "[brandes]")
{
	std::valarray<double> globalBC = { 1,2,2,1.5,1,3 };
	std::vector<std::shared_ptr<VertexInfo<int, double>>> verticesInfo(globalBC.size());
	for (int i = 0; i < verticesInfo.size(); ++i)
	{
		verticesInfo[i] = std::make_shared<VertexInfo<int, double>>(3);
	}
	verticesInfo[0]->setBorderSPLength(0, 1.0f);
	verticesInfo[0]->setBorderSPLength(1, 2.0f);
	verticesInfo[0]->setBorderSPLength(2, 3.0f);
	verticesInfo[0]->setBorderSPCount(0, 2);
	verticesInfo[0]->setBorderSPCount(1, 1);
	verticesInfo[0]->setBorderSPCount(2, 1);

	verticesInfo[1]->setBorderSPLength(0, 2.0f);
	verticesInfo[1]->setBorderSPLength(1, 1.0f);
	verticesInfo[1]->setBorderSPLength(2, 3.0f);
	verticesInfo[1]->setBorderSPCount(0, 2);
	verticesInfo[1]->setBorderSPCount(1, 2);
	verticesInfo[1]->setBorderSPCount(2, 1);

	verticesInfo[2]->setBorderSPLength(0, 2.0f);
	verticesInfo[2]->setBorderSPLength(1, 3.0f);
	verticesInfo[2]->setBorderSPLength(2, 4.0f);
	verticesInfo[2]->setBorderSPCount(0, 2);
	verticesInfo[2]->setBorderSPCount(1, 1);
	verticesInfo[2]->setBorderSPCount(2, 1);

	verticesInfo[3]->setBorderSPLength(0, 4.0f);
	verticesInfo[3]->setBorderSPLength(1, 3.0f);
	verticesInfo[3]->setBorderSPLength(2, 5.0f);
	verticesInfo[3]->setBorderSPCount(0, 2);
	verticesInfo[3]->setBorderSPCount(1, 2);
	verticesInfo[3]->setBorderSPCount(2, 1);

	verticesInfo[4]->setBorderSPLength(0, 5.0f);
	verticesInfo[4]->setBorderSPLength(1, 1.0f);
	verticesInfo[4]->setBorderSPLength(2, 3.0f);
	verticesInfo[4]->setBorderSPCount(0, 1);
	verticesInfo[4]->setBorderSPCount(1, 1);
	verticesInfo[4]->setBorderSPCount(2, 3);

	std::valarray<double> verticesClassCardinality(globalBC.size());

	std::set<int> vertices = { 0,1,2,3,4 };

	VertexInfoPivotSelector<int, double> ps;

	std::vector<int> pivots = 
		ps.selectPivots(globalBC, verticesInfo, verticesClassCardinality, vertices);

	REQUIRE(pivots.size() == 3);

	REQUIRE(std::find(pivots.begin(), pivots.end(), 0) != pivots.end());
	REQUIRE(std::find(pivots.begin(), pivots.end(), 3) != pivots.end());
	REQUIRE(std::find(pivots.begin(), pivots.end(), 4) != pivots.end());

	REQUIRE(verticesClassCardinality[0] == 2);
	REQUIRE(verticesClassCardinality[1] == 2);
	REQUIRE(verticesClassCardinality[2] == 2);
	REQUIRE(verticesClassCardinality[3] == 2);
	REQUIRE(verticesClassCardinality[4] == 1);
}