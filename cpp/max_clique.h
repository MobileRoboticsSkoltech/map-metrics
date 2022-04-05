#ifndef MAP_METRICS_MAX_CLIQUE_H
#define MAP_METRICS_MAX_CLIQUE_H

#include <vector>

#include "clustering.h"

namespace max_clique{
    struct MaxCliqueVisitor{
        MaxCliqueVisitor(clustering::ClusterMeans const & clusterizer);

        template<typename Clique, typename Graph>
        void clique(Clique const & c, Graph const & g){
            
        }


        clustering::ClusterMeans clusterizer_;
        size_t max_size_;
        std::vector<Eigen::Index> max_clique_idx_; 
    };

    void FindMaxClique(clustering::ClusterMeans const & clusterizer, double eps);
}

#endif // MAP_METRICS_MAX_CLIQUE_H
