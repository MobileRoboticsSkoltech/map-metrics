// Copyright (c) 2022, Skolkovo Institute of Science and Technology (Skoltech)
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//
//  Created on: May 20, 2022
//       Author: Arthur Saliou
//               arthur.salio@gmail.com
//
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
