#include <catch2/catch.hpp>

#include <DirectedWeightedGraph.h>

#include <exception>
#include <fstream>
#include <memory>

using namespace fastbc;

TEST_CASE("Directed weighted graph constructor/getters", "[fastbc]")
{
	std::ifstream dwgText("DWGtext.txt");
	if (!dwgText.is_open())
	{
		throw std::runtime_error("Unable to read test graph file.");
	}

	std::shared_ptr<IGraph<int, double>> graph;

	REQUIRE_NOTHROW(graph = std::make_shared<DirectedWeightedGraph<int, double>>(dwgText));

	REQUIRE(graph->vertices().size() == 9);
	REQUIRE(graph->edges() == 16);

	const auto& fs = graph->forwardStar(4);
	REQUIRE(fs.size() == 3);
	REQUIRE(fs.find(5)->second == 1);
	REQUIRE(fs.find(6)->second == 5);
	REQUIRE(fs.find(8)->second == 3);

	auto& bs = graph->backwardStar(4);
	REQUIRE(bs.size() == 3);
	REQUIRE(bs.find(0)->second == 7);
	REQUIRE(bs.find(2)->second == 4);
	REQUIRE(bs.find(3)->second == 3);

	REQUIRE(graph->edge(7, 5) == 2);
	REQUIRE(graph->edge(0, 1) == 4);
	REQUIRE(graph->edge(1, 0) == 0);
}