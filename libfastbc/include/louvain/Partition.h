#ifndef FASTBC_LOUVAIN_PARTITION_H
#define FASTBC_LOUVAIN_PARTITION_H

#include <louvain/LouvainGraph.h>
#include <algorithm>

namespace fastbc {
	namespace louvain {
		template <typename V, typename W>
		class Partition {
		public:
			LouvainGraph<V,W> g; // network to compute communities for
			int size; // nummber of nodes in the network and size of all std::vectors

			std::vector<double> neigh_weight;
			std::vector<unsigned int> neigh_pos;
			unsigned int neigh_last;
			std::vector<int> n2c; // community to which each node belongs

            //used to compute the modularity of the graph and the gains for each single node
			std::vector<double> woutc,      // sum of weights going out from node i (index) to nodes in his community
                                winc,       // sum of weights going in to node i (index) from nodes in his community
                                wout,       // weighted out degree of node i
                                win,        // weighted in degree of node i
                                woutctot,   // sum of all weights of edges exiting from nodes in community i (index)
                                winctot;    // sum of all weights of edges entering into nodes in community i (index)

			// number of pass for one level computation
			// if -1, compute as many pass as needed to increase modularity
			int nb_pass;

			// a new pass is computed if the last one has generated an increase 
			// greater than min_modularity
			// if 0. even a minor increase is enough to go for one more pass
			double min_modularity;

			Partition(LouvainGraph<V, W>& gc, double minm)  {
			    g = gc;
			    size = g.nb_nodes;

			    neigh_weight.resize(size,-1);
			    neigh_pos.resize(size);
			    neigh_last=0;

			    n2c.resize(size);
                woutc.resize(size);
                winc.resize(size);
                wout.resize(size);
			    win.resize(size);
			    woutctot.resize(size);
                winctot.resize(size);

			    for (int i=0 ; i<size ; i++) {
			        n2c[i] = i;
			        woutc[i] = winc[i] = g.weighted_selfloops(i);
                    wout[i] = woutctot[i] = g.weighted_out_degree(i);
                    win[i] = winctot[i] = g.weighted_in_degree(i);
			    }

			    nb_pass = -1;
			    min_modularity = minm;
			}

			// remove the node from its current community with which it has dnodecomm links
			inline void remove(int node);

			// insert the node in comm with which it shares dnodecomm links
			inline void insert(int node, int comm);

			// compute the gain of modularity if node where inserted in comm
			inline double modularity_gain(int node, int comm, double wic);

			// compute the set of neighboring communities of node
			// for each community, gives the number of links from node to comm
			void neigh_comm(unsigned int node);

			// compute the modularity of the current partition
			double modularity();

			// generates the binary graph of communities as computed by one_level
			LouvainGraph<V, W> partition2graph();

			// compute communities of the graph for one level
			// return true if some nodes have been moved
			bool one_level();
			bool one_level(std::vector<int> evaluation_order);

            void write_communities();
        };
    }
}

template<typename V, typename W>
inline void
fastbc::louvain::Partition<V, W>::remove(int node) {
  woutc[node] = 0;
  winc[node] = 0;
  woutctot[n2c[node]] -= g.weighted_out_degree(node);
  winctot[n2c[node]] -= g.weighted_in_degree(node);
  n2c[node]  = -1;
}

template <typename V, typename W>
inline void
fastbc::louvain::Partition<V, W>::insert(int node, int comm) {
  n2c[node]=comm;
  std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator > pin = g.out_neighbors(node);
  for (unsigned int i=0 ; i<g.nb_out_neighbors(node) ; i++) {
      if(n2c[(V)*(pin.first+i)] == comm)
        woutc[node] += (W)*(pin.second+i);
  }
  std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator > pout = g.in_neighbors(node);
  for (unsigned int i=0 ; i<g.nb_in_neighbors(node) ; i++) {
      if(n2c[(V)*(pout.first+i)] == comm)
        winc[node] += (W)*(pout.second+i);
  }
  woutctot[comm] += g.weighted_out_degree(node);
  winctot[comm] += g.weighted_out_degree(node);
}

template <typename V, typename W>
inline double
fastbc::louvain::Partition<V, W>::modularity_gain(int node, int comm, double wic) {
  double woutn  = wout[node];
  double winn   = win[node];
  double winc   = winctot[comm];
  double woutc  = woutctot[comm];
  double m      = (double) g.total_weight;

  /*std::cout << "wic   : " << wic << std::endl;
  std::cout << "woutn : " << woutn << std::endl;
  std::cout << "winn  : " << winn << std::endl;
  std::cout << "winc  : " << winc << std::endl;
  std::cout << "woutc : " << woutc << std::endl;
  std::cout << "m     : " << m << std::endl;
  std::cout << "gain  : " << (wic/m - (woutn/m)*(winc/m) - (winn/m)*(woutc/m)) << std::endl;*/
  
  return (wic/m - (woutn/m)*(winc/m) - (winn/m)*(woutc/m));
}

template<typename V, typename W>
double fastbc::louvain::Partition<V, W>::modularity() {
    double q    = 0.;
    double m = (double)g.total_weight;
    for (int i=0 ; i<size ; i++) {
        if (wout[i]>0){
            q += (double)woutc[i]/m - ((double)wout[i]/m)*((double)winctot[n2c[i]]/m);
        }
    }
    return q;
}

template<typename V, typename W>
void fastbc::louvain::Partition<V, W>::neigh_comm(unsigned int node) {
    for (unsigned int i=0 ; i<neigh_last ; i++)
        neigh_weight[neigh_pos[i]]=-1;
    neigh_last=0;

    std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator> pin     = g.in_neighbors(node);
    std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator> pout    = g.out_neighbors(node);

    unsigned int indeg  = g.nb_in_neighbors(node);
    unsigned int outdeg = g.nb_out_neighbors(node);

    neigh_pos[0]=n2c[node];
    neigh_weight[neigh_pos[0]]=0;
    neigh_last=1;

    for (unsigned int i=0 ; i<indeg ; i++) {
        unsigned int neigh          = *(pin.first+i);
        unsigned int neigh_comm     = n2c[neigh];
        double neigh_w = (g.inweights.size()==0)?1.:*(pin.second+i);
        
        if (neigh!=node) {
            if (neigh_weight[neigh_comm]==-1) {
                	neigh_weight[neigh_comm]=0.;
                	neigh_pos[neigh_last++]=neigh_comm;
            }
            neigh_weight[neigh_comm]+=neigh_w;
        }
    }

    for (unsigned int i=0 ; i<outdeg ; i++) {
        unsigned int neigh          = *(pout.first+i);
        if(neigh != node) {
            unsigned int neigh_comm     = n2c[neigh];
            double neigh_w = (g.outweights.size()==0)?1.:*(pout.second+i);
            
            if (neigh_weight[neigh_comm]==-1) {
                    neigh_weight[neigh_comm]=0.;
                    neigh_pos[neigh_last++]=neigh_comm;
            }
            neigh_weight[neigh_comm]+=neigh_w;
        }
    }
}


template<typename V, typename W>
fastbc::louvain::LouvainGraph<V, W> fastbc::louvain::Partition<V, W>::partition2graph() {
    // Renumber communities
    std::vector<int> renumber(size, -1);
    for (int node=0 ; node<size ; node++) {
        renumber[n2c[node]]++;
    }

    int final=0;
    for (int i=0 ; i<size ; i++)
        if (renumber[i]!=-1)
            renumber[i]=final++;

    // Compute communities
    std::vector<std::vector<int> > comm_nodes(final);
    for (int node=0 ; node<size ; node++) {
        comm_nodes[renumber[n2c[node]]].push_back(node);
    }

    // Compute weighted graph
    LouvainGraph<V, W> g2;
    g2.nb_nodes = comm_nodes.size();
    g2.total_weight = 0;
    g2.nb_links = 0;
    g2.indegrees.resize(comm_nodes.size());
    g2.outdegrees.resize(comm_nodes.size());

    int comm_deg = comm_nodes.size();
    for (int comm=0 ; comm<comm_deg ; comm++) {
        std::map<int,double> inm;
        std::map<int, double> outm;
        std::map<int,double>::iterator it;

        int comm_size = comm_nodes[comm].size();
        for (int node=0 ; node<comm_size ; node++) {
            std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator> pin;
            pin = g.in_neighbors(comm_nodes[comm][node]);
            int indeg = g.nb_in_neighbors(comm_nodes[comm][node]);
            for (int i=0 ; i<indeg ; i++) {
            	int neigh           = *(pin.first+i);
            	int neigh_comm      = renumber[n2c[neigh]];
            	double neigh_weight = (g.inweights.size()==0)?1.:*(pin.second+i);

            	it = inm.find(neigh_comm);
            	if (it==inm.end())
            	    inm.insert(std::make_pair(neigh_comm, neigh_weight));
            	else
            	    it->second+=neigh_weight;
            }

            std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator> pout;
            pout = g.out_neighbors(comm_nodes[comm][node]);
            int outdeg = g.nb_out_neighbors(comm_nodes[comm][node]);
            for (int i=0 ; i<outdeg ; i++) {
                int neigh           = *(pout.first+i);
                int neigh_comm      = renumber[n2c[neigh]];
                double neigh_weight = (g.outweights.size()==0)?1.:*(pout.second+i);

                it = outm.find(neigh_comm);
                if (it==outm.end())
                    outm.insert(std::make_pair(neigh_comm, neigh_weight));
                else
                    it->second+=neigh_weight;
            }
        }
        g2.indegrees[comm]=(comm==0)? inm.size() : g2.indegrees[comm-1]+inm.size();
        g2.outdegrees[comm]=(comm==0)? outm.size() : g2.outdegrees[comm-1]+outm.size();
        g2.nb_links+=inm.size();

        
        for (it = inm.begin() ; it!=inm.end() ; it++) {
            g2.total_weight    += it->second;
            g2.inlinks.push_back(it->first);
            g2.inweights.push_back(it->second);
        }

        for (it = outm.begin() ; it!=outm.end() ; it++) {
            g2.outlinks.push_back(it->first);
            g2.outweights.push_back(it->second);
        }
    }

    return g2;
}

/*
 * Maybe use deterministic seed.
 */
template<typename V, typename W>
bool fastbc::louvain::Partition<V, W>::one_level() {
    std::vector<int> random_order(size);
    for (int i=0 ; i<size ; i++)
        random_order[i]=i;
    std::random_shuffle(random_order.begin(), random_order.end());    
    return one_level(random_order);
}

template<typename V, typename W>
bool fastbc::louvain::Partition<V, W>::one_level(std::vector<int> random_order) {
    bool improvement=false ;
    int nb_moves;
    int nb_pass_done = 0;
    double new_mod     = modularity();
    double cur_mod     = new_mod;

    // repeat while 
    //     there is an improvement of modularity
    //     or there is an improvement of modularity greater than a given epsilon 
    //     or a predefined number of pass have been done
    do {
        cur_mod = new_mod;

        //cerr << "========== Iteration ==========    " << endl;
        //cerr << "modularity: " << cur_mod << endl << endl;

        nb_moves = 0;
        nb_pass_done++;

        // for each node: remove the node from its community and insert it in the best community
        for (int node_tmp=0 ; node_tmp<size ; node_tmp++) {
//            int node = node_tmp;
            int node = random_order[node_tmp];

            //cerr << "Node " << node << ": " << endl << endl;
            int node_comm         = n2c[node];

            // computation of all neighboring communities of current node
            neigh_comm(node);
            // remove node from its current community
            remove(node);

            // compute the nearest community for node
            // default choice for future insertion is the former community
            int best_comm        = node_comm;
            double best_nblinks  = 0.;
            double best_increase = 0.;
            for (unsigned int i=0 ; i<neigh_last ; i++) {
                //cerr << "Evaluating move to " << neigh_pos[i] << " community" << endl;
                double increase = modularity_gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]]);
                //cerr << "Gain: " << increase << endl;
                if (increase>best_increase) {
                    best_comm     = neigh_pos[i];
                    best_nblinks  = neigh_weight[neigh_pos[i]];
                    best_increase = increase;
                }
            }
            //cerr << "Best gain: " << best_increase << endl << endl;
            //cerr << "Inserting into " << best_comm << " community." << endl;
            // insert node in the nearest community

            insert(node, best_comm);
         
            //cerr << "New modularity: " << modularity() << endl << endl;

            if (best_comm!=node_comm)
                nb_moves++;
        }

        new_mod = modularity();
        //cerr << "Final iteration modularity: " << new_mod << endl;
        //cerr << "======== End Iteration ========    " << endl;
        if (nb_moves>0)
            improvement=true;
        
    } while (nb_moves>0 && new_mod-cur_mod>min_modularity);

    /*std::cout << "Communities: " << std::endl;
    write_communities();*/
    
    return improvement;
}

template<typename V, typename W>
void fastbc::louvain::Partition<V, W>::write_communities() {
  std::map<int, std::vector<int> > comms;
  for(int i=0; i<size; i++) {
    if(comms.count(n2c[i]) <= 0)
      comms[n2c[i]] = std::vector<int>();
    comms[n2c[i]].push_back(i);
  }

  std::map<int, std::vector<int> >::iterator it;

  for ( it = comms.begin(); it != comms.end(); it++ )
  {
      std::cout << it->first << " : " << std::endl;
      for(int i=0; i< it->second.size(); i++)
        std::cout << it->second[i] << " ";
      std::cout << std::endl;
  }
}

#endif