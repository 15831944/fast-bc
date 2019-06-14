#include "popl.hpp"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#define FASTBC_BRANDES_ENABLE_PIVOT_BORDER
#define FASTBC_BRANDES_CLUSTERED_IGNORE_UNCONNECTED

#include <DirectedWeightedGraph.h>
#include <brandes/ClusteredBrandesBC.h>
#include <brandes/DijkstraClusterEvaluator.h>
#include <brandes/DijkstraSSBrandesBC.h>
#include <brandes/ExactBrandesBC.h>
#include <brandes/KMeansPivotSelector.h>
#include <brandes/VertexInfoPivotSelector.h>
#include <kmeans/PlusPlusKMeans.h>
#include <louvain/LouvainGraphPartition.h>

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

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#endif

#define FASTBC_SPDLOG_FORMAT_DEBUG	"[%H:%M:%S.%f] %^[%=9l]%$ [%=7t] [%!]\n\t%v"
#define FASTBC_SPDLOG_FORMAT		"%^[%=9l]%$ %v"

int main(int argc, char **argv)
{
	/*
	 *	Program options 
	 */
	std::string edgeListPath, outBCPath, louvainSeed, loggerLevel;
	int louvainExecutors;
	double louvainPrecision, kFrac;
	bool exactBC;

	popl::OptionParser op("Usage: fastbc [ options ] <edge_list_path>");
	auto ls = op.add<popl::Value<std::string>, popl::Attribute::optional>(
		"s", "louvain-seeds",
		"Seeds to be used by each parallel louvain execution",
		"");
	ls->assign_to(&louvainSeed);
	auto le = op.add<popl::Value<int>, popl::Attribute::optional>(
		"e", "louvain-instances",
		"Number of parallel louvain instances",
		4);
	le->assign_to(&louvainExecutors);
	op.add<popl::Value<double>, popl::Attribute::optional>(
		"p", "louvain-precision",
		"Minimum precision value for louvain algorithm",
		0.01,
		&louvainPrecision);
	auto kf = op.add<popl::Value<double>, popl::Attribute::optional>(
		"k", "kfrac",
		"Topological classes aggregation factor (0-1). Enables 2-Clustered Brandes algorithm");
	kf->assign_to(&kFrac);
	op.add<popl::Switch, popl::Attribute::optional>(
		"", "exact",
		"Force exact betweenness computation (very long time)",
		&exactBC);
	op.add<popl::Value<std::string>, popl::Attribute::optional>(
		"o", "output",
		"Output file path",
		"bc.txt",
		&outBCPath);
	op.add<popl::Value<std::string>, popl::Attribute::optional>(
		"d", "debug",
		"Logger level (trace|debug|info|warning|error|critical|off)",
		"info",
		&loggerLevel
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
		std::cout << "Missing input file path" << "\n\n" << op.help();
		return -1;
	}
	else
	{
		edgeListPath = op.non_option_args().front();
	}

	// Setup logger
	spdlog::set_default_logger(spdlog::stdout_color_mt("fastbc"));
	auto log_level = spdlog::level::from_str(loggerLevel);
	if(log_level <= spdlog::level::debug)
	{
		spdlog::set_pattern(FASTBC_SPDLOG_FORMAT_DEBUG);
	}
	else
	{
		spdlog::set_pattern(FASTBC_SPDLOG_FORMAT);
	}
	spdlog::set_level(log_level);

	// Check bc output file
	std::ifstream outFileTest(outBCPath, std::ifstream::in);
	if (outFileTest.good())
	{
		SPDLOG_CRITICAL("File \"{}\" already existing", outBCPath);
		return -2;
	}
	outFileTest.close();

	// Initialize louvain seeds
	std::set<std::mt19937::result_type> seed;
	if (ls->is_set())
	{
		if (!le->is_set())
		{
			SPDLOG_CRITICAL("Louvain executors count must be set to allow executors seeds to be set.");
			return -1;
		}

		std::stringstream ss(louvainSeed);

		std::mt19937::result_type s;
		while (ss >> s)
		{
			if (!seed.insert(s).second)
			{
				SPDLOG_CRITICAL("Duplicate value in louvain seeds, each seed must be unique.");
				return -1;
			}

			if (ss.peek() == ',')
			{
				ss.ignore();
			}
		}

		if (seed.size() != louvainExecutors)
		{
			SPDLOG_CRITICAL("Louvain seeds count is different from louvain executors count.");
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

	// Check kfrac value range
	if (kf->is_set())
	{
		if (kFrac <= 0.0 || kFrac >= 1.0)
		{
			SPDLOG_CRITICAL("Kfrac value must be in range 0-1.");
			return -1;
		}
	}

	/*
	 *	Program options end
	 */



	/*
	 *	Program initialization
	 */
	// Open graph text file
	std::ifstream graphTextFile(edgeListPath);
	if (!graphTextFile.is_open())
	{
		SPDLOG_CRITICAL("There was an error opening given edge list file path.");
		return -1;
	}

	// Initialize graph object with loaded text file
	std::shared_ptr<fastbc::IGraph<FASTBC_V_TYPE, FASTBC_W_TYPE>> graph =
		std::make_shared<fastbc::DirectedWeightedGraph<FASTBC_V_TYPE, FASTBC_W_TYPE>>(graphTextFile);

	// Print some information about loaded graph
	SPDLOG_INFO("Loaded graph contains {} vertices and {} edges", graph->vertices().size(), graph->edges());

	std::shared_ptr<fastbc::brandes::IBrandesBC<FASTBC_V_TYPE, FASTBC_W_TYPE>> brandesBC;
	if(exactBC)
	{
		SPDLOG_INFO("Algorithm: exact Brandes' betweenness centrality");
		brandesBC = 
			std::make_shared<fastbc::brandes::ExactBrandesBC<FASTBC_V_TYPE, FASTBC_W_TYPE>>();
	}
	else
	{
		/* Louvain community detector */
		std::shared_ptr<fastbc::IGraphPartition<FASTBC_V_TYPE, FASTBC_W_TYPE>> louvainEvaluator =
			std::make_shared<fastbc::louvain::LouvainGraphPartition<FASTBC_V_TYPE, FASTBC_W_TYPE>>(
				seed, louvainPrecision);

		/* Brandes cluster evaluator */
		std::shared_ptr<fastbc::brandes::IClusterEvaluator<FASTBC_V_TYPE, FASTBC_W_TYPE>> clusterEvaluator =
			std::make_shared<fastbc::brandes::DijkstraClusterEvaluator<FASTBC_V_TYPE, FASTBC_W_TYPE>>();

		/* Cluster pivot selector */
		std::shared_ptr<fastbc::brandes::IPivotSelector<FASTBC_V_TYPE, FASTBC_W_TYPE>> pivotSelector;
		if (kf->is_set())
		{
			SPDLOG_INFO("Algorithm: 2-clustered Brandes' betweenness centrality");
			// Kmeans approximated pivot selector
			pivotSelector = 
				std::make_shared<fastbc::brandes::KMeansPivotSelector<FASTBC_V_TYPE, FASTBC_W_TYPE>>(
					std::shared_ptr<fastbc::brandes::IPivotSelector<FASTBC_V_TYPE, FASTBC_W_TYPE>>(
						new fastbc::brandes::VertexInfoPivotSelector<FASTBC_V_TYPE, FASTBC_W_TYPE>()),
					std::shared_ptr<fastbc::kmeans::IKMeans<FASTBC_V_TYPE, FASTBC_W_TYPE>>(
						new fastbc::kmeans::PlusPlusKMeans<FASTBC_V_TYPE, FASTBC_W_TYPE>()),
					kFrac);
		}
		else
		{
			SPDLOG_INFO("Algorithm: clustered Brandes' betweenness centrality");
			pivotSelector = 
				std::make_shared<fastbc::brandes::VertexInfoPivotSelector<FASTBC_V_TYPE, FASTBC_W_TYPE>>();
		}

		/* Single source Brandes */
		std::shared_ptr<fastbc::brandes::DijkstraSSBrandesBC<FASTBC_V_TYPE, FASTBC_W_TYPE>> singleSourceBC =
			std::make_shared<fastbc::brandes::DijkstraSSBrandesBC<FASTBC_V_TYPE, FASTBC_W_TYPE>>();

		/* Clustered Brandes Betweenness centrality calculator */
		brandesBC =
			std::make_shared<fastbc::brandes::ClusteredBrandeBC<FASTBC_V_TYPE, FASTBC_W_TYPE>>(
				louvainEvaluator, clusterEvaluator, singleSourceBC, pivotSelector);
	}
	

	/*
	 *	Program initialization end
	 */


	auto startTime = std::chrono::high_resolution_clock::now();

	std::vector<FASTBC_W_TYPE> bc = brandesBC->computeBC(graph);

	auto totalTime = std::chrono::high_resolution_clock::now() - startTime;
	auto milliTime = std::chrono::duration_cast<std::chrono::milliseconds>(totalTime).count();
	auto microTime = std::chrono::duration_cast<std::chrono::microseconds>(totalTime).count() - milliTime * 1000;

	SPDLOG_INFO("Total computation time: {}.{}ms", milliTime, microTime);

	/*
	 *	Save results
	 */
	std::ofstream outFile(outBCPath, std::ofstream::out);
	for (size_t i = 0; i < bc.size(); ++i)
	{
		if(bc[i] >= 0)
		{
			outFile << bc[i] << std::endl;
		}
		else
		{
			outFile << 0 << std::endl;
		}
		
	}

	SPDLOG_INFO("Results written to \"{}\"", outBCPath);

	return 0;
}
