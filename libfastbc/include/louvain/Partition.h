#ifndef FASTBC_LOUVAIN_PARTITION_H
#define FASTBC_LOUVAIN_PARTITION_H

#include <louvain/LouvainGraph.h>

namespace fastbc {
	namespace louvain {
		template <typename V, typename W>
		class Partition {
		public:

			std::vector<double> neigh_weight;
			std::vector<unsigned int> neigh_pos;
			unsigned int neigh_last;

			LouvainGraph<V,W> g; // network to compute communities for
			int size; // nummber of nodes in the network and size of all std::vectors
			std::vector<int> n2c; // community to which each node belongs
			std::vector<double> in,tot; // used to compute the modularity participation of each community

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
			    in.resize(size);
			    tot.resize(size);

			    for (int i=0 ; i<size ; i++) {
			        n2c[i] = i;
			        in[i]  = g.nb_selfloops(i);
			        tot[i] = g.weighted_degree(i);
			    }

			    nb_pass = -1;
			    min_modularity = minm;
			}

			// remove the node from its current community with which it has dnodecomm links
			inline void remove(int node, int comm, double dnodecomm);

			// insert the node in comm with which it shares dnodecomm links
			inline void insert(int node, int comm, double dnodecomm);

			// compute the gain of modularity if node where inserted in comm
			// given that node has dnodecomm links to comm.  The formula is:
			// [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
			// [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
			// where In(comm)    = number of half-links strictly inside comm
			//       Tot(comm)   = number of half-links inside or outside comm (sum(degrees))
			//       d(node,com) = number of links from node to comm
			//       deg(node)   = node degree
			//       m           = number of links
			inline double modularity_gain(int node, int comm, double dnodecomm, double w_degree);

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
		};

	}
}

template<typename V, typename W>
inline void
fastbc::louvain::Partition<V, W>::remove(int node, int comm, double dnodecomm) {
  tot[comm] -= g.weighted_degree(node);
  in[comm]  -= 2*dnodecomm + g.nb_selfloops(node);
  n2c[node]  = -1;
}

template <typename V, typename W>
inline void
fastbc::louvain::Partition<V, W>::insert(int node, int comm, double dnodecomm) {
  tot[comm] += g.weighted_degree(node);
  in[comm]  += 2*dnodecomm + g.nb_selfloops(node);
  n2c[node]=comm;
}

template <typename V, typename W>
inline double
fastbc::louvain::Partition<V, W>::modularity_gain(int node, int comm, double dnodecomm, double w_degree) {
  double totc = (double)tot[comm];
  double degc = (double)w_degree;
  double m2   = (double)g.total_weight;
  double dnc  = (double)dnodecomm;

  /*cerr << "tot      : ";
  for(int i=0; i<tot.size(); i++) cerr << tot[i] << " ";
  cerr << endl;
  cerr << "DiC      : " << dnc << endl;
  cerr << "Di       : " << degc << endl;
  cerr << "EtotC    : " << totc << endl;
  cerr << "2m       : " << m2 << endl;*/
  
  return (dnc - totc*degc/m2);
}

template<typename V, typename W>
double fastbc::louvain::Partition<V, W>::modularity() {
    double q    = 0.;
    double m2 = (double)g.total_weight;
    for (int i=0 ; i<size ; i++) {
        if (tot[i]>0){
            q += (double)in[i]/m2 - ((double)tot[i]/m2)*((double)tot[i]/m2);
        }
    }
    return q;
}

template<typename V, typename W>
void fastbc::louvain::Partition<V, W>::neigh_comm(unsigned int node) {
    for (unsigned int i=0 ; i<neigh_last ; i++)
        neigh_weight[neigh_pos[i]]=-1;
    neigh_last=0;

    std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator> p = g.neighbors(node);

    unsigned int deg = g.nb_neighbors(node);

    neigh_pos[0]=n2c[node];
    neigh_weight[neigh_pos[0]]=0;
    neigh_last=1;

    for (unsigned int i=0 ; i<deg ; i++) {
        unsigned int neigh                = *(p.first+i);
        unsigned int neigh_comm     = n2c[neigh];
        double neigh_w = (g.weights.size()==0)?1.:*(p.second+i);
        
        if (neigh!=node) {
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
    g2.degrees.resize(comm_nodes.size());

    int comm_deg = comm_nodes.size();
    for (int comm=0 ; comm<comm_deg ; comm++) {
        std::map<int,double> m;
        std::map<int,double>::iterator it;

        int comm_size = comm_nodes[comm].size();
        for (int node=0 ; node<comm_size ; node++) {
            std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator> p;
            p = g.neighbors(comm_nodes[comm][node]);
            int deg = g.nb_neighbors(comm_nodes[comm][node]);
            for (int i=0 ; i<deg ; i++) {
            	int neigh           = *(p.first+i);
            	int neigh_comm      = renumber[n2c[neigh]];
            	double neigh_weight = (g.weights.size()==0)?1.:*(p.second+i);

            	it = m.find(neigh_comm);
            	if (it==m.end())
            	    m.insert(std::make_pair(neigh_comm, neigh_weight));
            	else
            	    it->second+=neigh_weight;
            }
        }
        g2.degrees[comm]=(comm==0)? m.size() : g2.degrees[comm-1]+m.size();
        g2.nb_links+=m.size();

        
        for (it = m.begin() ; it!=m.end() ; it++) {
            g2.total_weight    += it->second;
            g2.links.push_back(it->first);
            g2.weights.push_back(it->second);
        }
    }

    return g2;
}

template<typename V, typename W>
bool fastbc::louvain::Partition<V, W>::one_level() {
    std::vector<int> random_order(size);
        for (int i=0 ; i<size ; i++)
                random_order[i]=i;
        /*for (int i=0 ; i<size-1 ; i++) {
                int rand_pos = rand()%(size-i)+i;
                int tmp            = random_order[i];
                random_order[i] = random_order[rand_pos];
                random_order[rand_pos] = tmp;
        }*/
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
            double w_degree = g.weighted_degree(node);

            // computation of all neighboring communities of current node
            neigh_comm(node);
            // remove node from its current community
            remove(node, node_comm, neigh_weight[node_comm]);

            // compute the nearest community for node
            // default choice for future insertion is the former community
            int best_comm                = node_comm;
            double best_nblinks    = 0.;
            double best_increase = 0.;
            for (unsigned int i=0 ; i<neigh_last ; i++) {
                //cerr << "Evaluating move to " << neigh_pos[i] << " community" << endl;
                double increase = modularity_gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]], w_degree);
                //cerr << "Gain: " << increase << endl;
                if (increase>best_increase) {
                    best_comm         = neigh_pos[i];
                    best_nblinks    = neigh_weight[neigh_pos[i]];
                    best_increase = increase;
                }
            }
            //cerr << "Best gain: " << best_increase << endl << endl;
            //cerr << "Inserting into " << best_comm << " community." << endl;
            // insert node in the nearest community
            insert(node, best_comm, best_nblinks);
         
            //cerr << "New modularity: " << modularity() << endl << endl;

            if (best_comm!=node_comm)
                nb_moves++;
        }

        double total_tot=0;
        double total_in=0;
        for (unsigned int i=0 ; i<tot.size() ;i++) {
            total_tot+=tot[i];
            total_in+=in[i];
        }

        new_mod = modularity();
        //cerr << "Final iteration modularity: " << new_mod << endl;
        //cerr << "======== End Iteration ========    " << endl;
        if (nb_moves>0)
            improvement=true;
        
    } while (nb_moves>0 && new_mod-cur_mod>min_modularity);

    return improvement;
}

#endif