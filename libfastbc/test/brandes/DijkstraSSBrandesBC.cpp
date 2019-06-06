#include <catch2/catch.hpp>

#include <brandes/DijkstraSSBrandesBC.h>

#include <DirectedWeightedGraph.h>
#include <SubGraph.h>
#include <fstream>
#include <memory>

using namespace fastbc::brandes;

TEST_CASE("Single source Brandes BC", "[brandes]")
{
	std::ifstream dwgText("DWGtext.txt");
	if (!dwgText.is_open())
	{
		throw std::exception("Unable to read test graph file.");
	}

	std::shared_ptr<fastbc::IGraph<int, float>> fullGraph =
		std::make_shared<fastbc::DirectedWeightedGraph<int, float>>(dwgText);

	std::shared_ptr<ISSBrandesBC<int, float>> ssBC =
		std::make_shared<DijkstraSSBrandesBC<int, float>>();

	std::valarray<float> globalBC = ssBC->singleSourceBrandes(0, fullGraph);

	REQUIRE(globalBC.size() == fullGraph->vertices().size());
}