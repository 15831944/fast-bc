#include <catch2/catch.hpp>

#include <SubGraph.h>

#include <DirectedWeightedGraph.h>
#include <exception>
#include <fstream>
#include <memory>

using namespace fastbc;

TEST_CASE("SubGraph contructor and getters", "[fastbc]")
{
	std::ifstream dwgText("DWGtext.txt");
	if (!dwgText.is_open())
	{
		throw std::exception("Unable to read test graph file.");
	}

	std::shared_ptr<IGraph<int, double>> graph;

	REQUIRE_NOTHROW(graph = std::make_shared<DirectedWeightedGraph<int, double>>(dwgText));

	std::shared_ptr<ISubGraph<int, double>> subGraph;

	REQUIRE_THROWS(subGraph = std::make_shared<SubGraph<int, double>>(std::set<int>({ 0, 1, 2, 3, 6 }), graph));

	REQUIRE_NOTHROW(subGraph = std::make_shared<SubGraph<int, double>>(std::set<int>({ 0, 1, 2, 3, 4 }), graph));

	REQUIRE(subGraph->vertices().size() == 5);
	REQUIRE(subGraph->edges() == 7);

	REQUIRE(subGraph->forwardStar(4).size() == 0);
	REQUIRE(subGraph->backwardStar(4).size() == 3);

	const auto& fs = subGraph->forwardStar(3);
	REQUIRE(fs.size() == 1);
	REQUIRE(fs.find(4)->second == 3);

	REQUIRE(subGraph->borders().size() == 2);
	REQUIRE(subGraph->isBorder(3));
	REQUIRE(subGraph->isBorder(4));

	REQUIRE(subGraph->referenceGraph() == graph);
}