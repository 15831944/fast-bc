#ifndef FASTBC_IDEGREEGRAPH_H
#define FASTBC_IDEGREEGRAPH_H

#include <IGraph.h>

namespace fastbc {
    template<typename V, typename W>
    class IDegreeGraph : public IGraph<V, W>
    {
    public:

        virtual void addEdge(V from, V to, W weight) = 0;

        virtual void initVertices() = 0;

        virtual W totalWeight() const = 0;

        virtual W inWeightedDegree(V v) const = 0;

        virtual W outWeightedDegree(V v) const = 0;
    };    
}


#endif