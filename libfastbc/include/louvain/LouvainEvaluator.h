#ifndef FASTBC_LOUVAIN_LOUVAINEVALUATOR_H
#define FASTBC_LOUVAIN_LOUVAINEVALUATOR_H

#include <louvain/ILouvainEvaluator.h>
#include <louvain/LouvainGraph.h>
#include <louvain/Partition.h>
#include <louvain/Community.h>
#include <ctime>
#include <omp.h>

namespace fastbc {
	namespace louvain {
			
		template<typename V, typename W>
		class LouvainEvaluator : public ILouvainEvaluator<V, W>
		{
	 
		private:
			typedef std::vector<std::shared_ptr<ICommunity<V,W>>> Result;
			typedef std::shared_ptr<IGraph<V,W>> Graph;
			bool verbose = true;	
			double precision = 0.01;	
			int parallelism = 4;

			void
			display_time(const char *str) {
				time_t rawtime;
				time ( &rawtime );
				std::cout << str << ": " << ctime (&rawtime);
			}

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
			build_result(const std::vector<int>& n2c, const Graph g) {
				Result r;

				int max = 0;
				for(int i=0; i<n2c.size(); i++)
					if(n2c[i] > max)
						max = n2c[i];

				r.resize(max + 1);
				for(int i=0; i<r.size(); i++)
					r[i] = std::make_shared<Community<V, W>>(g);

				for(int i=0; i<n2c.size(); i++)
					r[n2c[i]]->add(i);

				return r;
			}


		public:
			LouvainEvaluator() {}		
				
			Result evaluateGraph(Graph graph) override
			{
			    time_t time_begin, time_end;
			    time(&time_begin);
			    if (verbose)
			        display_time("Begin");

			    LouvainGraph<V, W> g(graph);
			    Partition<V, W> p(g, precision);

			    bool improvement=true;
			    double mod=p.modularity(), new_mod;
			    int level=0;

			    do {
			        if (verbose) {
			            std::cout << "level " << level << ":\n";
			            display_time("    start computation");
			            std::cout << "    network size: " 
				     << p.g.nb_nodes << " nodes, " 
				     << p.g.nb_links << " links, "
				     << p.g.total_weight << " weight." << std::endl;
			        }

			        improvement = p.one_level();
			        new_mod = p.modularity();
			        g = p.partition2graph();
			        renumber_communities(n2c, p.n2c);
			        p = Partition(g, precision);

			        if (verbose)
			            std::cout << "  modularity increased from " << mod << " to " << new_mod << std::endl;

			        mod=new_mod;
			        if (verbose)
			            display_time("  end computation");
			       level++;
			    } while(improvement);

			    time(&time_end);
			    if (verbose) {
			        display_time("End");
			        std::cout << "Total duration: " << (time_end-time_begin) << " sec." << std::endl; 
			    }
			    std::cout << new_mod << std::endl;
			    std::cout << std::endl;

				return build_result(n2c, graph);
			}
		};
	}
}


#endif
