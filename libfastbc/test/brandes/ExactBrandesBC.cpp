#include <catch2/catch.hpp>

#include <brandes/ExactBrandesBC.h>

#include <DirectedWeightedGraph.h>
#include <SubGraph.h>
#include <fstream>

using namespace fastbc::brandes;

TEST_CASE("Exact Brandes' BC computation test", "[brandes]")
{
    std::ifstream dwgText("DWGtext.txt");
	if (!dwgText.is_open())
	{
		throw std::runtime_error("Unable to read test graph file.");
	}

    std::shared_ptr<fastbc::IGraph<int, float>> fullGraph = 
		std::make_shared<fastbc::DirectedWeightedGraph<int, float>>(dwgText);

	std::shared_ptr<fastbc::ISubGraph<int, float>> subGraph =
		std::make_shared<fastbc::SubGraph<int, float>>(std::vector<int>({ 0,1,2,3,4 }), fullGraph);

    std::shared_ptr<IBrandesBC<int, float>> exactBrandesBC = 
        std::make_shared<ExactBrandesBC<int, float>>();

    std::vector<float> graphBC = exactBrandesBC->computeBC(subGraph);

    // Check betweenness centrality values
	REQUIRE(graphBC[0] == 2.0f);
	REQUIRE(graphBC[1] == 2.0f);
	REQUIRE(graphBC[2] == 0.5f);
	REQUIRE(graphBC[3] == 1.0f);
	REQUIRE(graphBC[4] == 0.0f);
}