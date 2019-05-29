#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <utility>

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

	parse_args(argc, argv);

	//**************************************************************************************
	// READ GRAPH

    return 0;
}