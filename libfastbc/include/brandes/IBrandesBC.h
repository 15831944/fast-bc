#ifndef FASTBC_BRANDES_IBRANDESBC_H
#define FASTBC_BRANDES_IBRANDESBC_H

#include <IGraph.h>

#include <memory>
#include <vector>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class IBrandesBC
		{
		public:

			virtual std::vector<W> computeBC(const std::shared_ptr<const IGraph<V, W>> graph) = 0;
		};

	}
}

#endif