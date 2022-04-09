#include "max_clique.h"
#include "clustering.h"

#include <memory>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>

#include <iostream>

namespace max_clique{
    // TODO: Use pointer instead of copying clusterizer since it can be big.
    MaxCliqueVisitor::MaxCliqueVisitor(clustering::ClusterMeans const & clusterizer,
                    std::shared_ptr<Eigen::Index> max_covered_points, 
                    std::shared_ptr<std::vector<Eigen::Index>> max_clique_idx)
                : clusterizer_(clusterizer), max_covered_points_(max_covered_points), 
                  max_clique_idx_(max_clique_idx){
    }

    std::vector<Eigen::Index> FindMaxClique(clustering::ClusterMeans const & clusterizer, double eps){
        typedef boost::adjacency_matrix<boost::undirectedS> Graph; 
        
        Graph adj_mx(clusterizer.getClusterNumber());
        for (Eigen::Index i = 0; i < clusterizer.getClusterNumber(); ++i){
            for (Eigen::Index j = 0; j < i; ++j){
                if (std::abs(clusterizer.getMeans().col(i)
                        .dot(clusterizer.getMeans().col(j))) < eps){
                    boost::add_edge(i, j, adj_mx);
                    boost::add_edge(j, i, adj_mx);            
                }
            }
        }

        auto max_covered_points = std::make_shared<Eigen::Index>(0);
        auto max_clique_idx = std::make_shared<std::vector<Eigen::Index>>();

        MaxCliqueVisitor vis(clusterizer, max_covered_points, max_clique_idx);
        boost::bron_kerbosch_all_cliques(adj_mx, vis);

        return *max_clique_idx;
    }
}
