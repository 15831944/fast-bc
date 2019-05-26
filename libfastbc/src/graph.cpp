#include <graph.h>

using namespace fastbc;

Graph::Graph(char *filename) {
	std::ifstream finput;
	finput.open(filename,std::fstream::in);	
  	if (finput.is_open() != true) {
    		std::cerr << "The file " << filename << " does not exist" << std::endl;
    		exit(EXIT_FAILURE);
  	}

  	unsigned long long nb_links = 0ULL;

  	while (!finput.eof()) {
    		unsigned int src, dest;
    		long double weight = 1.0L;

    		finput >> src >> dest >> weight;
    
    		if (finput) {
      			if (links.size()<= std::max(src,dest)+1) {
        			links.resize(std::max(src,dest)+1);
      			}
      
      			links[src].push_back(std::make_pair(dest,weight));
      			if (src!=dest)
        			links[dest].push_back(std::make_pair(src,weight));

      			nb_links += 1ULL;
    		}
  	}

  	finput.close();
}