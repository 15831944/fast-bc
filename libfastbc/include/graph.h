#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

namespace fastbc 
{

	class Graph 
	{
	public:
		std::vector<std::vector<std::pair<int, long double> > > links;
	  
		Graph (char *filename);
	};

}


#endif // GRAPH_H
