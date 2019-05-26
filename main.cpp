#include <graph.h>

int main(int argc, char *argv)
{
	fastbc::Graph g("../graph.txt");
	std::cout << g.links.size() << std::endl;
	return 1;
}