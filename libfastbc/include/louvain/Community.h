#ifndef FASTBC_LOUVAIN_COMMUNITY_H
#define FASTBC_LOUVAIN_COMMUNITY_H

#include <louvain/ICommunity.h>
#include <set>

namespace fastbc {

	namespace louvain {

		template<typename V, typename W>
		class Community : public ICommunity<V,W>
		{
		public:
			Community(std::shared_ptr<const IGraph<V, W>> g) : g(g) {}

			Community(std::set<V> u, std::shared_ptr<const IGraph<V, W>> g) : u(u), g(g) {}

			const std::set<V>& all() const override 
			{
				return u;
			}

			bool contains(V vertex) const override
			{
				return u.find(vertex) != u.end();
			}

			void add(V vertex) override
			{
				u.insert(vertex);
			}

			void remove(V vertex) override
			{
				u.erase(vertex);
			}

			std::shared_ptr<const IGraph<V, W>> referenceGraph() const override
			{
				return g;
			}

			int size() const override 
			{
				return u.size();
			}
		private:
			std::set<V> u;
			std::shared_ptr<const IGraph<V, W>> g;
		};
	}
}
#endif
