#include <DirectedWeightedGraph.h>
#include <SubGraph.h>
#include <louvain/LouvainEvaluator.h>
#include <louvain/Community.h>

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
	
	fastbc::louvain::ILouvainEvaluator<int, float>* l = new fastbc::louvain::LouvainEvaluator<int, float>();
	
	std::vector<std::shared_ptr<fastbc::louvain::ICommunity<int, float>>> p = l->evaluateGraph(graph);
	/*for(int i=0; i<p.comms.size(); i++) {
		std::cout << "Community: " << i << std::endl;
		std::vector<int> vec = p.comms[i]->all();
		for(int j=0; j<vec.size(); j++)
			std::cout << vec[j] << " ";
		std::cout << std::endl;
	}*/

	// TODO: Print computation stats

	// TODO: Print result to file

	return 0;
}
