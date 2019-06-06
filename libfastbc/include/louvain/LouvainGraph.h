#ifndef FASTBC_LOUVAIN_LOUVAINGRAPH_H
#define FASTBC_LOUVAIN_LOUVAINGRAPH_H

namespace fastbc {
	namespace louvain {

		template<typename V, typename W>
		class LouvainGraph {
		 public:
			typedef std::shared_ptr<const IDegreeGraph<V,W>> Graph;
			
			unsigned int nb_nodes;
			unsigned long nb_links;
			W total_weight;    

			std::vector<unsigned long> indegrees;
			std::vector<V> inlinks;
			std::vector<W> inweights;

            std::vector<unsigned long> outdegrees;
            std::vector<V> outlinks;
            std::vector<W> outweights;

			LouvainGraph() {}
			LouvainGraph(Graph graph);

			// return the number of neighbors (degree) of the node
			inline unsigned int nb_in_neighbors(V node);
            inline unsigned int nb_out_neighbors(V node);

			// return the number of self loops of the node
			inline W weighted_selfloops(V node);

			// return the weighted degree of the node
			inline W weighted_in_degree(V node);
            inline W weighted_out_degree(V node);

			// return pointers to the first neighbor and first weight of the node
			inline std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator > in_neighbors(V node);
            inline std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator > out_neighbors(V node);
		};

	}
}

template<typename V, typename W>
fastbc::louvain::LouvainGraph<V, W>::LouvainGraph(Graph graph) {
    nb_nodes = graph->vertices().size();
    nb_links = graph->edges();

    indegrees.resize(nb_nodes);
    inlinks.resize(nb_links);
    inweights.resize(nb_links);

    outdegrees.resize(nb_nodes);
    outlinks.resize(nb_links);
    outweights.resize(nb_links);

    unsigned long intot = 0, outtot = 0;
    total_weight = 0;
    for(int i=0; i<nb_nodes; i++) {
        for(auto &[v, w]: graph->backwardStar(i)) {
            inlinks[intot] = v;
            inweights[intot] = w;
            intot += 1;
            total_weight += w;
        }
        indegrees[i] = intot;

        for(auto &[v, w]: graph->forwardStar(i)) {
            outlinks[outtot] = v;
            outweights[outtot] = w;
            outtot += 1;
        }
        outdegrees[i] = outtot;
    }
}

template<typename V, typename W>
inline unsigned int
fastbc::louvain::LouvainGraph<V, W>::nb_in_neighbors(V node) {
    if (node==0)
        return indegrees[0];
    else
        return indegrees[node]-indegrees[node-1];
}

template<typename V, typename W>
inline unsigned int
fastbc::louvain::LouvainGraph<V, W>::nb_out_neighbors(V node) {
    if (node==0)
        return outdegrees[0];
    else
        return outdegrees[node]-outdegrees[node-1];
}

template<typename V, typename W>
inline W
fastbc::louvain::LouvainGraph<V, W>::weighted_selfloops(V node) {
    std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator > p = out_neighbors(node);
    for (unsigned int i=0 ; i<nb_out_neighbors(node) ; i++) {
        if (*(p.first+i)==node) {
            if (outweights.size()!=0)
				return (W)*(p.second+i);
            else 
				return 1.;
        }
    }
    return 0.;
}

template<typename V, typename W>
inline W
fastbc::louvain::LouvainGraph<V, W>::weighted_in_degree(V node) {
    if (inweights.size()==0)
        return (W)nb_in_neighbors(node);
    else {
        std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator > p = in_neighbors(node);
        W res = 0;
        for (unsigned int i=0 ; i<nb_in_neighbors(node) ; i++) {
            res += (W)*(p.second+i);
        }
        return res;
    }
}

template<typename V, typename W>
inline W
fastbc::louvain::LouvainGraph<V, W>::weighted_out_degree(V node) {
    if (outweights.size()==0)
        return (W)nb_out_neighbors(node);
    else {
        std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator > p = out_neighbors(node);
        W res = 0;
        for (unsigned int i=0 ; i<nb_out_neighbors(node) ; i++) {
            res += (W)*(p.second+i);
        }
        return res;
    }
}

template<typename V, typename W>
inline std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator >
fastbc::louvain::LouvainGraph<V, W>::in_neighbors(V node) {
    if (node==0)
        return std::make_pair(inlinks.begin(), inweights.begin());
    else if (inweights.size()!=0)
        return std::make_pair(inlinks.begin()+indegrees[node-1], inweights.begin()+indegrees[node-1]);
    else
        return std::make_pair(inlinks.begin()+indegrees[node-1], inweights.begin());
}

template<typename V, typename W>
inline std::pair<typename std::vector<V>::iterator, typename std::vector<W>::iterator >
fastbc::louvain::LouvainGraph<V, W>::out_neighbors(V node) {
    if (node==0)
        return std::make_pair(outlinks.begin(), outweights.begin());
    else if (outweights.size()!=0)
        return std::make_pair(outlinks.begin()+outdegrees[node-1], outweights.begin()+outdegrees[node-1]);
    else
        return std::make_pair(outlinks.begin()+outdegrees[node-1], outweights.begin());
}


#endif