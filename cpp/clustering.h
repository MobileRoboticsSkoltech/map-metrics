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
#ifndef MAP_METRICS_CLUSTERING_H
#define MAP_METRICS_CLUSTERING_H

#include <cstdint>
#include <vector>
#include <cilantro/core/data_containers.hpp>
#include <alglib/ap.h>

namespace clustering{

    class ClusterMeans{
     public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        explicit ClusterMeans(alglib::integer_1d_array const & labels, Eigen::Index cluster_number);

        void filterClusters(Eigen::Ref<const cilantro::VectorSet3d> const points, int32_t min_clust_size);
        
        const cilantro::VectorSet3d & getMeans() const;

        const Eigen::VectorXi & getIdx() const;

        const Eigen::VectorXi & getLabels() const;

        Eigen::Index getClusterNumber() const;

     private:
        Eigen::VectorXi labels_;
        cilantro::VectorSet3d cluster_means_;
        Eigen::VectorXi cluster_idx_;
        Eigen::Index cluster_number_;
    };

    ClusterMeans ClusterizeAHC(Eigen::Ref<const cilantro::VectorSet3d> const points, 
                                   double distance_treshold);

}

#endif // MAP_METRICS_CLUSTERING_H