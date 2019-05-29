#ifndef FASTBC_IGRAPH_H
#define FASTBC_IGRAPH_H

#include <map>
#include <vector>

namespace fastbc {

    /**
     * @brief Generic grpah object class
     * 
     * @tparam V Type for vertex index number
     * @tparam W Type for edge weight value
     */
    template<typename V, typename W>
    class IGraph
    {
    public:

        /**
         * @brief Get weight of given src->dest edge
         * 
         * @details When src->dest edge is present its weight will be returned, else zero 
         *          will be the return value. 
         *          In case of undirected graph the calls edge(src,dest) == edge(dest,src).
         * 
         * @param src Source node for required edge
         * @param dest Destination node for required edge
         * @return W Weight of required edge
         */
        virtual W edge(V src, V dest) const = 0;

        /**
         * @brief Get forward star vertex/weight for given src vertex
         * 
         * @param src Vertex index 
         * @return const std::map<V, W>& Dest/edge weight map of all outgoing edges from src vertex
         */
        virtual const std::map<V, W>& forwardStar(V src) const = 0;

        virtual V nodes() const = 0;

        virtual V edges() const = 0;
    };    

}


#endif