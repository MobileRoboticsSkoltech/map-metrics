#ifndef MAP_METRICS_MAX_CLIQUE_H
#define MAP_METRICS_MAX_CLIQUE_H

#include <vector>
#include <memory>

#include "clustering.h"
#include "boost/graph/adjacency_matrix.hpp"
#include <iostream>

namespace max_clique{
    class MaxCliqueVisitor{
     public:
        MaxCliqueVisitor(clustering::ClusterMeans const & clusterizer,
                        std::shared_ptr<Eigen::Index> max_covered_points, 
                        std::shared_ptr<std::vector<Eigen::Index>> max_clique_idx);

        template<typename Clique, typename Graph>
        void clique(Clique const & c, Graph const & g){
            if (c.size() == 2) return; 
            
            std::vector<Eigen::Index> clique_idx;
            clique_idx.reserve(c.size());

            Eigen::Index number_of_points = 0;

            typename Clique::const_iterator i, end = c.end();
            for (i = c.begin(); i != end; ++i){
                number_of_points += (clusterizer_.getLabels().array() == clusterizer_.getIdx()[*i]).count();
                clique_idx.push_back(*i);
            }

            if (number_of_points > *max_covered_points_){
                *max_covered_points_ = number_of_points;
                *max_clique_idx_ = clique_idx;
            }
        }

     private:
        clustering::ClusterMeans clusterizer_;
        std::shared_ptr<Eigen::Index> max_covered_points_;
        std::shared_ptr<std::vector<Eigen::Index>> max_clique_idx_; 
    };

    std::vector<Eigen::Index> FindMaxClique(clustering::ClusterMeans const & clusterizer, double eps);
}

#endif // MAP_METRICS_MAX_CLIQUE_H
