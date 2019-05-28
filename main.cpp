#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/property_map/shared_array_property_map.hpp>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <utility>

// TODO choose the right data structures 
// https://www.boost.org/doc/libs/1_55_0/libs/graph/doc/using_adjacency_list.html

using namespace boost;

char *filename 		= NULL;
char *filename_o 	= NULL;
bool verbose 		= false;

void usage(char *prog_name, const char *more) 
{
    std::cerr << more;
    std::cerr << "usage: " << prog_name << " input_file [-o output_file] [-v] [-h]" << std::endl << std::endl;
    std::cerr << "input_file: file containing the graph." << std::endl;
    std::cerr << "-o file\tfile to write output betweennes centrality.";
    std::cerr << "-v\tverbose mode." << std::endl;
    std::cerr << "-h\tshow this usage message." << std::endl;
    exit(0);
}

void parse_args(int argc, char **argv) {
    if (argc<2)
        usage(argv[0], "Bad arguments number\n");

    for (int i = 1; i < argc; i++) {
        if(argv[i][0] == '-') {
            switch(argv[i][1]) {
            case 'o':
                filename_o = argv[i+1];
                i++;
                break;
            case 'v':
            	verbose=true;
            	break;
            case 'h':
            	usage(argv[0], "\n");
            default:
	            usage(argv[0], "Unknown option\n");
            }
        } else {
            if (filename==NULL)
                filename = argv[i];
            else
                usage(argv[0], "More than one filename\n");
        }
    }
}

int main(int argc, char **argv)
{
	srand ( time(NULL) );
	typedef adjacency_list<vecS, listS, directedS,
		property<vertex_index_t, int>,
		property<edge_capacity_t, double> 
	> Graph;

	parse_args(argc, argv);

	//**************************************************************************************
	// READ GRAPH

	Graph g;
    auto ensure_node = [&g, known=std::map<int, Graph::vertex_descriptor>{}](int id) mutable {
        if (auto it = known.find(id); it == known.end())
            return known.emplace(id, add_vertex(id, g)).first->second;
        else
            return it->second;
    };

    std::ifstream ifs(filename);

    for (std::string s; getline(ifs, s);) {
        int from,to;
        double weight;
        if (std::istringstream(s) >> from >> to >> weight) {
            add_edge(ensure_node(from), ensure_node(to), weight, g);
        }
        else break;
    }
    ifs.close();

    if (ifs.eof()) {
        std::cout << "Graph parsed with " << num_edges(g) << " edges and " << num_vertices(g) << " vertices\n";
    } else {
        ifs.clear();
        ifs.close();
        std::cout << "Parse error\n";
        exit(EXIT_FAILURE);
    }

    if (ifs) {
        std::cout << "Remaining unparsed input: '" << ifs.rdbuf() << "'\n";
        ifs.close();
        exit(EXIT_FAILURE);
    } else {
   		ifs.close();
    }


    //**************************************************************************************
	// PERFORM BRANDES

	typedef property_map< Graph, vertex_index_t>::type VertexIndexMap;
	VertexIndexMap v_index = get(vertex_index, g);
	boost::shared_array_property_map<double, VertexIndexMap >
		centrality_map(num_vertices(g), v_index);

    brandes_betweenness_centrality(g, centrality_map);

    BGL_FORALL_VERTICES_T(v, g, Graph) {
    	std::cout << v_index[v] << " " << get(centrality_map, v) << std::endl;
	}

    return 0;
}