#ifndef FASTBC_LOUVAIN_LOUVAINGRAPHPARTITION_H
#define FASTBC_LOUVAIN_LOUVAINGRAPHPARTITION_H

#include <IGraphPartition.h>
#include <louvain/LouvainGraph.h>
#include <louvain/Partition.h>

#include <memory>
#include <random>
#include <spdlog/spdlog.h>

namespace fastbc {
	namespace louvain {
			
		template<typename V, typename W>
		class LouvainGraphPartition : public IGraphPartition<V, W>
		{
	 
		private:
			typedef std::vector<std::vector<V>> Result;
			typedef std::shared_ptr<const IDegreeGraph<V,W>> Graph;

			double _precision;	
			int _parallelism;
			std::vector<std::mt19937> _seed;

			void
			renumber_communities(std::vector<V>& comms, const std::vector<V>& n2c) {
				    std::vector<int> renumber(n2c.size(), -1);
				    for (int node=0 ; node<n2c.size() ; node++) {
				        renumber[n2c[node]]++;
				    }

				    int final=0;
				    for (int i=0 ; i<n2c.size() ; i++)
				        if (renumber[i]!=-1)
				            renumber[i]=final++;

				    for(int i=0; i<comms.size(); i++) {
				    	int old_community = comms[i];
				    	int new_community = n2c[old_community];
				    	int new_renumbered_community = renumber[new_community];

				    	comms[i] = new_renumbered_community;
				    }
			}

			Result
			build_result(const std::vector<int>& n2c, Graph g) {
				Result r;

				int max = 0;
				for(size_t i=0; i<n2c.size(); i++)
				{
					if(n2c[i] > max)
					{
						max = n2c[i];
					}
				}

				r.resize(max + 1);

				for(size_t i=0; i<n2c.size(); i++)
				{
					r[n2c[i]].push_back(i);
				}

				return r;
			}


		public:
			LouvainGraphPartition(
				const std::set<std::mt19937::result_type>& seeds, 
				double precision = 0.01)
				: _parallelism(seeds.size()), _precision(precision)
			{
				for (auto& seed : seeds)
				{
					_seed.push_back(std::mt19937(seed));
				}
			}		
				
			Result partitionGraph(Graph graph) override
			{
			    LouvainGraph<V, W> g(graph);
			    std::vector<Partition<V, W> > p(_parallelism, Partition<V, W>(g, _precision));
			    std::vector<V> n2c(g.nb_nodes);
			    for(int i=0; i<g.nb_nodes; i++) n2c[i] = i;
			    std::vector<bool> improvements(_parallelism, true);
			    std::vector<double> modularities(_parallelism);
			    int best_i = 0;

			    bool improvement;
			    double mod=p[best_i].modularity(), new_mod;
			    int level=0;

			    do {
					SPDLOG_DEBUG("Level: {}\n\tNetwork size: {} vertices, {} edges, {} weight",
						level, g.nb_nodes, g.nb_links, g.total_weight);

			        #pragma omp parallel for
			        for(int i=0; i<_parallelism; i++) {
			        	improvements[i] = p[i].one_level(_seed[i]);
			        	modularities[i] = p[i].modularity();
			        }

			        int best_i = 0;
			        for(int i=1; i<_parallelism; i++)
			        	if(modularities[i] > modularities[best_i])
			        		best_i = i;

			        improvement = improvements[best_i];
			        new_mod = p[best_i].modularity();
			        g = p[best_i].partition2graph();
			        renumber_communities(n2c, p[best_i].n2c);
			        for(int i=0; i<_parallelism; i++)
			        	p[i] = Partition<V, W> (g, _precision);

					SPDLOG_DEBUG("Modularity increased from {} to {}", mod, new_mod);

			        mod = new_mod;
					
					level++;
			    } while(improvement);

				SPDLOG_DEBUG("Final modularity {}", new_mod);

				return build_result(n2c, graph);
			}
		};
	}
}


#endif
