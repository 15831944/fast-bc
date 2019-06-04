#ifndef FASTBC_BRANDES_IPIVOTSELECTOR_H
#define FASTBC_BRANDES_IPIVOTSELECTOR_H

#include "VertexInfo.h"

#include <memory>
#include <set>
#include <vector>

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class IPivotSelector
		{
		public:

			virtual std::vector<V> selectPivots(
				const std::vector<std::shared_ptr<VertexInfo<V, W>>>& verticesList) = 0;

			virtual std::vector<V> selectPivots(
				const std::vector<std::shared_ptr<VertexInfo<V, W>>>& verticesList,
				const std::set<V>& vertices) = 0;
		};

	}
}

#endif // !FASTBC_BRANDES_IPIVOTSELECTOR_H
