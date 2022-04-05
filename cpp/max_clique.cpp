#include "max_clique.h"
#include "clustering.h"

#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>

namespace max_clique{
    void FindMaxClique(clustering::ClusterMeans const & clusterizer, double eps){
        typedef boost::adjacency_matrix<boost::undirectedS> Graph; 
        
        Graph adj_mx(clusterizer.getClusterNumber());
        for (Eigen::Index i = 0; i < clusterizer.getClusterNumber(); ++i){
            for (Eigen::Index j = 0; j < i; ++j){
                // TODO: Remove excessive copying of clusterMeans
                if (std::abs(clusterizer.getMeans().col(i)
                        .dot(clusterizer.getMeans().col(j))) < eps){
                    boost::add_edge(i, j, adj_mx);
                    boost::add_edge(j, i, adj_mx);            
                }
            }
        }

        // MaxCliqueVisitor vis(clusterizer);

    }
}
