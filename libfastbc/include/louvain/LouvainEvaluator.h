#ifndef FASTBC_LOUVAIN_LOUVAINEVALUATOR_H
#define FASTBC_LOUVAIN_LOUVAINEVALUATOR_H

#include <louvain/ILouvainEvaluator.h>
#include <louvain/Community.h>
#include <louvain/Partition.h>
#include <time.h>

namespace fastbc {
	namespace louvain {
			
		template<typename V, typename W>
		class LouvainEvaluator : public ILouvainEvaluator<V, W>
		{
	 
		private:
			typedef std::vector<std::shared_ptr<ICommunity<V,W>>> Result;
			typedef std::shared_ptr<IGraph<V,W>> Graph;
		
			bool verbose = false;

			Partition<V, W> singleton_partition(Graph g)
			{
				Partition<V, W> p;
				p.num_comms = 0;
				const std::vector<V>& v = g->verticesList();
				for (auto it = v.begin(); it != v.end(); ++it) {
					std::shared_ptr<ICommunity<V,W>> comm = 
						std::make_shared<Community<V,W>>(g);

					comm->add(*it);
					p.comms.push_back(comm);
					p.wd_in_for_comm.push_back(g->inWeightedDegree(*it));
					p.wd_out_for_comm.push_back(g->outWeightedDegree(*it));
					p.v2c.push_back(p.comms.size() - 1);
					p.num_comms += 1;
				}
				return p;
			}

			bool move_nodes(Graph g, Partition<V, W>& p)
			{
				double h_old, h_new = modularity(g, p), precision = 0.01;
				int n_moves;
				bool improvement = false;
				const std::vector<V>& v = g->verticesList();

				do {
					std::cout << "========== Iteration ==========  " << std::endl;
					if(verbose) {
    					std::cout << "modularity: " << h_new << std::endl << std::endl;
					}
					n_moves = 0;
					h_old = h_new;

					//TODO here random order
					for(int i=0; i<v.size(); i++) {
						if(verbose) {
							std::cout << "Node " << i << ": " << std::endl << std::endl;
						}
						int node_comm = p.v2c[v[i]];

						std::unordered_set<int> comms = neigh_comms(g, p, v[i]);

						remove_from_comm(g, p, v[i], node_comm);
						
						int best_comm = node_comm;
						double best_gain = 0.;
						for (auto const & c: comms) {
							double gain = delta_modularity(g, p, v[i], c);
							if(verbose) {
								std::cout << "Evaluating move to " << c << " community" << std::endl;
								std::cout << "Gain: " << gain << std::endl << std::endl;
							}
							if(gain > best_gain) {
								best_gain = gain;
								best_comm = c;
							}
					    }
					    if(verbose) {
					    	std::cout << "Best gain: " << best_gain << std::endl << std::endl;
      						std::cout << "Inserting into " << best_comm << " community." << std::endl;
      					}
					    add_to_comm(g, p, v[i], best_comm);

					    if(verbose) {
					    	std::cout << "New modularity: " << modularity(g, p) << std::endl << std::endl;
						} else {
							std::cout << "Evaluated node " << v[i] << std::endl;
						}

					    if(best_comm != node_comm)
					    	n_moves++;
					}
					h_new = modularity(g, p);

					if(verbose) {
						std::cout << "Final iteration modularity: " << h_new << std::endl;
    				}
    				std::cout << "======== End Iteration ========  " << std::endl;
					if (n_moves>0)
      					improvement=true;				
				}while(n_moves > 0 && h_new-h_old>precision);
				
				return improvement;
			}

			void remove_from_comm(Graph g, Partition<V, W>&p, V v, int c) {
				//Update old community
				p.comms[c]->remove(v);
				if(p.comms[c]->size() == 0)
					p.num_comms = p.num_comms - 1;
				p.wd_in_for_comm[c] -= g->inWeightedDegree(v);
				p.wd_out_for_comm[c] -= g->outWeightedDegree(v);
				p.v2c[v] = -1;
			}

			void add_to_comm(Graph g, Partition<V, W>& p, V v, int c) {
				p.v2c[v] = c;
				p.comms[c]->add(v);
				p.wd_in_for_comm[c] += g->inWeightedDegree(v);
				p.wd_out_for_comm[c] += g->outWeightedDegree(v);
			}

			double modularity(Graph g, Partition<V, W> p) {
				double m = g->totalWeight();
				double s = 0;
				for(int i=0; i<p.v2c.size(); i++) {
					double woutc = compute_weighted_degree_to_community(g, p, i, p.v2c[i]);
					double wout = g->outWeightedDegree(i);
					double wintotc = p.wd_in_for_comm[p.v2c[i]];
					s += woutc - wout*wintotc/m;
				}
				double q = s/m;
				return q;
			}

			inline double compute_weighted_degree_to_community(Graph g, Partition<V, W> p, V v, int c) {
				double wd = 0;
				for ( const auto &[outv, outw]: g->forwardStar(v) )
					if(p.v2c[outv] == c)
						wd += outw;
				return wd;
			}

			inline double compute_weighted_degree_from_community(Graph g, Partition<V, W> p, V v, int c) {
				double wd = 0;
				for ( const auto &[inv, inw]: g->backwardStar(v) )
					if(p.v2c[inv] == c)
						wd += inw;
				return wd;
			}

			inline double compute_weighted_degree_in_community(Graph g, Partition<V, W> p, V v, int c) {
				double wd = 0;
				wd += compute_weighted_degree_from_community(g, p, v, c);
				wd += compute_weighted_degree_to_community(g, p, v, c);
				return wd;
			}


			double delta_modularity(Graph g, Partition<V, W> p, V v, int c) {
				double m = g->totalWeight();
				double wc = compute_weighted_degree_in_community(g, p, v, c);
				double wout = g->outWeightedDegree(v);
				double win = g->inWeightedDegree(v);
				double einc = p.wd_in_for_comm[c];
				double eoutc = p.wd_out_for_comm[c];
				if(verbose) {				
					std::cout << "WiC      : " << wc << std::endl;
					std::cout << "Wouti    : " << wout << std::endl;
					std::cout << "Wini     : " << win << std::endl;
					std::cout << "EinC     : " << einc << std::endl;
					std::cout << "EoutC    : " << eoutc << std::endl;
					std::cout << "m        : " << m << std::endl;
				}
				double d = wc - (wout*einc + win*eoutc)/m;
				return d/m;
			}

			std::unordered_set<int> neigh_comms(Graph g, Partition<V, W> p, V v) {
				std::unordered_set<int> u;

				u.insert(p.v2c[v]);

				std::map<V, W> forwardStar = g->forwardStar(v);
				std::map<V, W> backwardStar = g->backwardStar(v);

				for ( const auto &[key, value]: forwardStar )
					u.insert(p.v2c[key]);
				for ( const auto &[key, value]: backwardStar )
					u.insert(p.v2c[key]);

				return u;
			}

			Graph aggregate_graph(Graph g, Partition<V, W> p)
			{
				std::map<V, std::map<V, W>> edges;

				Graph new_graph = std::make_shared<fastbc::DirectedWeightedGraph<V, W>>();

				//Iterate all edges by community
				//for each community
				for(int i=0; i<p.comms.size(); i++) {
					edges[i] = std::map<V, W>();
					//for each vertex in the community
					std::vector<V> v = p.comms[i]->all();
					for(int j=0; j<v.size(); j++) {
						//get every outgoing edge;
						std::map<V, W> fws = g->forwardStar(v[j]);
						//foreach outgoing edge
						for ( const auto &[out_v, out_w]: fws ) {
							//update the corresponding weight
							if(edges[i].count(out_v) <= 0) {
								edges[i][out_v] = out_w;
							} else {
								edges[i][out_v] += out_w;
							}
						}
					}
				}

				for ( const auto &[from_v, to_vw]: edges )
					for ( const auto &[to_v, w]: to_vw )
						new_graph->addEdge(from_v, to_v, w);

				new_graph->initVertices();

				return new_graph;
			}

			Result init_result(Graph g) {
				Result result;
				const std::vector<V>& v = g->verticesList();
				for (auto it = v.begin(); it != v.end(); ++it) {
					std::shared_ptr<ICommunity<V,W>> comm = 
						std::make_shared<Community<V,W>>(g);
					comm->add(*it);
					result.push_back(comm);
				}
				return result;
			}

			Result update_result(Graph g, Partition<V, W> p, Result old_res) {
				Result res;
				for(int i=0; i<p.comms.size(); i++) {
					if(p.comms[i]->size() > 0) {
						std::shared_ptr<ICommunity<V,W>> comm = std::make_shared<Community<V,W>>(g);
						std::vector<V> v = p.comms[i]->all();
						for(int j=0; j<v.size(); j++) {
							std::vector<V> old_v = old_res[v[j]]->all();
							for(int k=0; k<old_v.size(); k++)
								comm->add(old_v[k]);
						}
						res.push_back(comm);
					}
				}
				return res;
			}

			void print_result(Result res) {
				std::cout << std::endl << "========== RESULT ==========" << std::endl << std::endl;
				for(int i=0; i<res.size(); i++) {
					std::cout << "Community " << i << ":" << std::endl;
					std::vector<V> v = res[i]->all();
					for(int j=0; j<v.size(); j++)
						std::cout << v[j] << " ";
					std::cout << std::endl << std::endl;
				}
			}

		public:
			LouvainEvaluator() {}		
				
			Result evaluateGraph(Graph graph) override
			{
				bool done = true;
				Graph g = graph;
				Partition<V, W> p = singleton_partition(g);
				Result res = init_result(g);

				double mod = modularity(g, p);
				do {
					done = move_nodes(g, p);
					mod = modularity(g, p);
					g = aggregate_graph(g, p);
					res = update_result(g, p, res);
					p = singleton_partition(g);
				} while(!done);

				if(verbose) {
					print_result(res);
				}

				std::cout << "Modularity:";
				std::cout << mod << std::endl;
				return res;
			}
		};
	}
}


#endif
