#ifndef FASTBC_LOUVAIN_COMMUNITY_H
#define FASTBC_LOUVAIN_COMMUNITY_H

#include <louvain/ICommunity.h>
#include <unordered_set>

namespace fastbc {

	namespace louvain {

		template<typename V, typename W>
		class Community : public ICommunity<V,W>
		{
		public:
			Community(std::shared_ptr<const IGraph<V, W>> g) : g(g)
			{
				u = std::unordered_set<V>();
			}

			Community(std::unordered_set<V> u, std::shared_ptr<const IGraph<V, W>> g) : u(u), g(g) {}

			std::vector<V> all() const override 
			{
				return std::vector(u.begin(), u.end());
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
			std::unordered_set<V> u;
			std::shared_ptr<const IGraph<V, W>> g;
		};
	}
}
#endif
