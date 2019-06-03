#ifndef FASTBC_LOUVAIN_LOUVAINGRAPH_H
#define FASTBC_LOUVAIN_LOUVAINGRAPH_H

namespace fastbc {
	namespace louvain {

		template<typename V, typename W>
		class LouvainGraph {
		 public:
			typedef std::shared_ptr<IGraph<V,W>> Graph;
			
			unsigned int nb_nodes;
			unsigned long nb_links;
			W total_weight;    

			std::vector<unsigned long> degrees;
			std::vector<V> links;
			std::vector<W> weights;

			LouvainGraph() {}
			LouvainGraph(Graph graph);

			// return the number of neighbors (degree) of the node
			inline unsigned int nb_neighbors(V node);

			// return the number of self loops of the node
			inline W nb_selfloops(V node);

			// return the weighted degree of the node
			inline W weighted_degree(V node);

			// return pointers to the first neighbor and first weight of the node
			inline std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator > neighbors(V node);
		};

	}
}

template<typename V, typename W>
fastbc::louvain::LouvainGraph<V, W>::LouvainGraph(std::shared_ptr<IGraph<V,W>> graph) {
    nb_nodes = graph->vertices();
    nb_links = graph->edges();

    degrees.resize(nb_nodes);
    links.resize(nb_links);
    weights.resize(nb_links);
    unsigned long tot = 0;
    total_weight = 0;
    for(int i=0; i<nb_nodes; i++) {
        for(auto &[v, w]: graph->forwardStar(i)) {
            links[tot] = v;
            weights[tot] = w;
            tot += 1;
            total_weight += w;
        }
        degrees[i] = tot;
    }
}

template<typename V, typename W>
inline unsigned int
fastbc::louvain::LouvainGraph<V, W>::nb_neighbors(V node) {
    if (node==0)
        return degrees[0];
    else
        return degrees[node]-degrees[node-1];
}

template<typename V, typename W>
inline W
fastbc::louvain::LouvainGraph<V, W>::nb_selfloops(V node) {
    std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator > p = neighbors(node);
    for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
        if (*(p.first+i)==node) {
            if (weights.size()!=0)
				return (W)*(p.second+i);
            else 
				return 1.;
        }
    }
    return 0.;
}

template<typename V, typename W>
inline W
fastbc::louvain::LouvainGraph<V, W>::weighted_degree(V node) {
    if (weights.size()==0)
        return (W)nb_neighbors(node);
    else {
        std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator > p = neighbors(node);
        W res = 0;
        for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
            res += (W)*(p.second+i);
        }
        return res;
    }
}

template<typename V, typename W>
inline std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator >
fastbc::louvain::LouvainGraph<V, W>::neighbors(V node) {
    if (node==0)
        return std::make_pair(links.begin(), weights.begin());
    else if (weights.size()!=0)
        return std::make_pair(links.begin()+degrees[node-1], weights.begin()+degrees[node-1]);
    else
        return std::make_pair(links.begin()+degrees[node-1], weights.begin());
}


#endif