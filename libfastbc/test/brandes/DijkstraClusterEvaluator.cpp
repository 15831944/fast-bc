#include <catch2/catch.hpp>

#include <brandes/DijkstraClusterEvaluator.h>

#include <DirectedWeightedGraph.h>
#include <SubGraph.h>
#include <exception>
#include <fstream>
#include <memory>
#include <valarray>
#include <vector>

using namespace fastbc::brandes;

TEST_CASE("Dijkstra cluster evaluation test", "[brandes]")
{
	std::ifstream dwgText("DWGtext.txt");
	if (!dwgText.is_open())
	{
		throw std::exception("Unable to read test graph file.");
	}

	std::shared_ptr<fastbc::IGraph<int, float>> fullGraph = 
		std::make_shared<fastbc::DirectedWeightedGraph<int, float>>(dwgText);

	std::shared_ptr<fastbc::ISubGraph<int, float>> subGraph =
		std::make_shared<fastbc::SubGraph<int, float>>(std::set<int>({ 0,1,2,3,4 }), fullGraph);

	std::shared_ptr<IClusterEvaluator<int, float>> ce = 
		std::make_shared<DijkstraClusterEvaluator<int, float>>();

	std::valarray<float> globalBC((float)0.0, fullGraph->vertices().size());
	std::vector<std::shared_ptr<VertexInfo<int, float>>> globalVertexInfo(fullGraph->vertices().size(), nullptr);

	ce->evaluateCluster(globalBC, globalVertexInfo, subGraph);

	// Check betweenness centrality values
	REQUIRE(globalBC[0] == 2.0f);
	REQUIRE(globalBC[1] == 2.0f);
	REQUIRE(globalBC[2] == 0.5f);
	REQUIRE(globalBC[3] == 1.0f);
	REQUIRE(globalBC[4] == 0.0f);

	for (int i = 5; i < globalBC.size(); ++i)
	{
		REQUIRE(globalBC[i] == 0.0f);
	}

	// Check vertices information

	REQUIRE(globalVertexInfo[0]->getBorderSPCount(0) == 1);
	REQUIRE(globalVertexInfo[0]->getBorderSPLength(0) == 5.0f);
	REQUIRE(globalVertexInfo[0]->getBorderSPCount(1) == 2);
	REQUIRE(globalVertexInfo[0]->getBorderSPLength(1) == 7.0f);

	REQUIRE(globalVertexInfo[1]->getBorderSPCount(0) == 1);
	REQUIRE(globalVertexInfo[1]->getBorderSPLength(0) == 1.0f);
	REQUIRE(globalVertexInfo[1]->getBorderSPCount(1) == 1);
	REQUIRE(globalVertexInfo[1]->getBorderSPLength(1) == 4.0f);

	REQUIRE(globalVertexInfo[2]->getBorderSPCount(0) == 1);
	REQUIRE(globalVertexInfo[2]->getBorderSPLength(0) == 8.0f);
	REQUIRE(globalVertexInfo[2]->getBorderSPCount(1) == 1);
	REQUIRE(globalVertexInfo[2]->getBorderSPLength(1) == 4.0f);

	REQUIRE(globalVertexInfo[3]->getBorderSPCount(0) == 1);
	REQUIRE(globalVertexInfo[3]->getBorderSPLength(0) == 0.0f);
	REQUIRE(globalVertexInfo[3]->getBorderSPCount(1) == 1);
	REQUIRE(globalVertexInfo[3]->getBorderSPLength(1) == 3.0f);

	REQUIRE(globalVertexInfo[4]->getBorderSPCount(0) == 0);
	REQUIRE(globalVertexInfo[4]->getBorderSPLength(0) == std::numeric_limits<float>::max());
	REQUIRE(globalVertexInfo[4]->getBorderSPCount(1) == 1);
	REQUIRE(globalVertexInfo[4]->getBorderSPLength(1) == 0.0f);
}