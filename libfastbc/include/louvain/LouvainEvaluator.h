#ifndef FASTBC_LOUVAIN_LOUVAINEVALUATOR_H
#define FASTBC_LOUVAIN_LOUVAINEVALUATOR_H

#include <louvain/ILouvainEvaluator.h>
#include <louvain/LouvainGraph.h>
#include <louvain/Partition.h>
#include <ctime>

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

			void
			display_time(const char *str) {
				time_t rawtime;
				time ( &rawtime );
				std::cout << str << ": " << ctime (&rawtime);
			}


		public:
			LouvainEvaluator() {}		
				
			Result evaluateGraph(Graph graph) override
			{
				Result r;
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
				return r;
			}
		};
	}
}


#endif
