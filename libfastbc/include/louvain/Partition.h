#ifndef FASTBC_LOUVAIN_PARTITION_H
#define FASTBC_LOUVAIN_PARTITION_H

#include <louvain/ICommunity.h>

namespace fastbc {
	namespace louvain {

		template <typename V, typename W>
		class Partition {
		public:
			std::vector<std::shared_ptr<ICommunity<V,W>>> comms;
			std::vector<int> v2c;
			std::vector<int> wd_in_for_comm;
			std::vector<int> wd_out_for_comm;
			int num_comms;
		};

	}
}

#endif