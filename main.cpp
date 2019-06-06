#include "popl.hpp"

#include <DirectedWeightedGraph.h>
#include <brandes/DijkstraClusterEvaluator.h>
#include <brandes/VertexInfoPivotSelector.h>
#include <brandes/DijkstraSSBrandesBC.h>
#include <brandes/ClusteredBrandesBC.h>
#include <louvain/LouvainEvaluator.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

#ifndef FASTBC_V_TYPE
#define FASTBC_V_TYPE int
#endif // !FASTBC_V_TYPE

#ifndef FASTBC_W_TYPE
#define FASTBC_W_TYPE double
#endif // !FASTBC_W_TYPE



int main(int argc, char **argv)
{
	std::string edgeListPath, outBCPath, louvainSeed;
	int louvainExecutors;
	double louvainPrecision;
	bool traceOn;

	popl::OptionParser op("Usage: fastbc [ options ] <edge_list_path>");
	auto ls = op.add<popl::Value<std::string>, popl::Attribute::optional>(
		"s", "louvain-seeds",
		"Seeds to be used by each parallel louvain execution",
		"");
	ls->assign_to(&louvainSeed);
	auto le = op.add<popl::Value<int>, popl::Attribute::optional>(
		"e", "louvain-executors",
		"Number of parallel louvain executor",
		4);
	le->assign_to(&louvainExecutors);
	op.add<popl::Value<double>, popl::Attribute::optional>(
		"p", "louvain-precision",
		"Minimum precision value for louvain algorithm",
		0.01,
		&louvainPrecision);
	op.add<popl::Value<std::string>, popl::Attribute::optional>(
		"o", "output",
		"Output file path",
		"bc.txt",
		&outBCPath);
	op.add<popl::Switch, popl::Attribute::optional>(
		"d", "debug",
		"Print centroids at each iteration step",
		&traceOn
		);

	try {
		op.parse(argc, argv);
	}
	catch (popl::invalid_option& e)
	{
		std::cout << e.what() << "\n\n" << op.help();
		return -1;
	}

	// Check if input file has been given
	if (op.non_option_args().size() != 1)
	{
		std::cout << "Missing edge list file path." << std::endl << std::endl << op.help();
		return -1;
	}
	else
	{
		edgeListPath = op.non_option_args().front();
	}

	// Check bc output file
	std::ifstream outFileTest(outBCPath, std::ifstream::in);
	if (outFileTest.good())
	{
		std::cout << "File " << outBCPath << " already existing." << std::endl;
		return -2;
	}
	outFileTest.close();

	// Initialize louvain seeds
	std::set<std::mt19937::result_type> seed;
	if (ls->is_set())
	{
		if (!le->is_set())
		{
			std::cout << "Louvain executors count must be set to allow executors seeds to be set." << std::endl;
			return -1;
		}

		std::stringstream ss(louvainSeed);

		std::mt19937::result_type s;
		while (ss >> s)
		{
			if (!seed.insert(s).second)
			{
				std::cout << "Duplicate value in louvain seeds, each seed must be unique." << std::endl;
				return -1;
			}

			if (ss.peek() == ',')
			{
				ss.ignore();
			}
		}

		if (seed.size() != louvainExecutors)
		{
			std::cout << "Louvain seeds count is different from louvain executors count." << std::endl;
			return -1;
		}
	}
	else
	{
		for (int i = 0; i < louvainExecutors; ++i)
		{
			seed.insert(std::chrono::high_resolution_clock::now().time_since_epoch().count());
		}
	}


	// Open graph text file
	std::ifstream graphTextFile(edgeListPath);
	if (!graphTextFile.is_open())
	{
		std::cout << "There was an error opening given edge list file path." << std::endl;
		return -1;
	}

	// Initialize graph object with loaded text file
	std::shared_ptr<fastbc::IGraph<FASTBC_V_TYPE, FASTBC_W_TYPE>> graph =
		std::make_shared<fastbc::DirectedWeightedGraph<FASTBC_V_TYPE, FASTBC_W_TYPE>>(graphTextFile);

	// Print some information about loaded graph
	std::cout << "Loaded graph contains " << graph->vertices().size() << " nodes and " 
		<< graph->edges() << " edges." << std::endl;


	std::shared_ptr<fastbc::louvain::ILouvainEvaluator<FASTBC_V_TYPE, FASTBC_W_TYPE>> louvainEvaluator =
		std::make_shared<fastbc::louvain::LouvainEvaluator<FASTBC_V_TYPE, FASTBC_W_TYPE>>(
			seed, louvainPrecision, traceOn);

	std::shared_ptr<fastbc::brandes::IClusterEvaluator<FASTBC_V_TYPE, FASTBC_W_TYPE>> clusterEvaluator =
		std::make_shared<fastbc::brandes::DijkstraClusterEvaluator<FASTBC_V_TYPE, FASTBC_W_TYPE>>();

	std::shared_ptr<fastbc::brandes::IPivotSelector<FASTBC_V_TYPE, FASTBC_W_TYPE>> pivotSelector =
		std::make_shared<fastbc::brandes::VertexInfoPivotSelector<FASTBC_V_TYPE, FASTBC_W_TYPE>>();

	std::shared_ptr<fastbc::brandes::DijkstraSSBrandesBC<FASTBC_V_TYPE, FASTBC_W_TYPE>> singleSourceBC =
		std::make_shared<fastbc::brandes::DijkstraSSBrandesBC<FASTBC_V_TYPE, FASTBC_W_TYPE>>();


	std::shared_ptr<fastbc::brandes::IBrandesBC<FASTBC_V_TYPE, FASTBC_W_TYPE>> brandesBC =
		std::make_shared<fastbc::brandes::ClusteredBrandeBC<FASTBC_V_TYPE, FASTBC_W_TYPE>>(
			louvainEvaluator, clusterEvaluator, singleSourceBC, pivotSelector);


	auto startTime = std::chrono::high_resolution_clock::now();

	std::valarray<FASTBC_W_TYPE> bc = brandesBC->computeBC(graph);

	auto totalTime = std::chrono::high_resolution_clock::now() - startTime;
	auto milliTime = std::chrono::duration_cast<std::chrono::milliseconds>(totalTime).count();
	auto microTime = std::chrono::duration_cast<std::chrono::microseconds>(totalTime).count() - milliTime * 1000;

	std::cout << "Total computation time: " << milliTime << "." << microTime << "ms\n";

	std::ofstream outFile(outBCPath, std::ofstream::out);
	for (size_t i = 0; i < bc.size(); ++i)
	{
		outFile << bc[i] << std::endl;
	}

	return 0;
}
