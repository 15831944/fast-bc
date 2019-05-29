#include <DirectedWeightedGraph.h>
#include <SubGraph.h>

#include <fstream>
#include <iostream>

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		return -1;
	}

	// TODO: Program options

	// Open graph text file
	std::ifstream graphTextFile(argv[1]);
	if (!graphTextFile.is_open())
	{
		std::cout << "There was an error opening given graph text file." << std::endl;
		return -1;
	}

	// Initialize graph object with loaded text file
	std::shared_ptr<fastbc::IGraph<int, float>> graph = 
		std::make_shared<fastbc::DirectedWeightedGraph<int, float>>(graphTextFile);

	// Print some information about loaded graph
	std::cout << "Loaded graph contains " << graph->vertices() << " nodes and " 
		<< graph->edges() << " edges." << std::endl;

	// TODO: Graph processing

	// TODO: Print computation stats

	// TODO: Print result to file

	return 0;
}