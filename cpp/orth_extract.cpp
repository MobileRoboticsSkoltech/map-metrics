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
#include "orth_extract.h"
#include "metrics_utils.h"
#include "clustering.h"
#include "max_clique.h"

#include <cilantro/core/kd_tree.hpp>
#include <alglib/dataanalysis.h>

#include <iostream>

namespace orth_extract{
    using PointCloud = cilantro::PointCloud3d;
    using KDTree = cilantro::KDTree3d<>;

    std::vector<cilantro::VectorSet3d> ExtractOrthogonalSubset(const PointCloud& pivot_pc, config::CustomConfig config, double eps){
        PointCloud cut_pc = EstimateNormals(pivot_pc, config.knn_rad, config.max_nn);
        auto normals = cut_pc.normals;

        clustering::ClusterMeans clusterizer = clustering::ClusterizeAHC(normals, 1e-1);
        clusterizer.filterClusters(normals, config.min_clust_size);

        std::vector<Eigen::Index> max_clique = max_clique::FindMaxClique(clusterizer, 1e-1);

        std::vector<cilantro::VectorSet3d> orthogonal_subset;
        orthogonal_subset.reserve(max_clique.size());

        for (Eigen::Index idx : max_clique){
            auto idx_mask = (clusterizer.getLabels().array() == clusterizer.getIdx()[idx]);
            cilantro::VectorSet3d new_points(3, idx_mask.count());
            Eigen::Index new_points_idx = 0;
            for (Eigen::Index i = 0; i < idx_mask.size(); ++i){
                if (idx_mask[i]) new_points.col(new_points_idx++) = cut_pc.points.col(i);
            }
            orthogonal_subset.push_back(new_points);
        }

        return orthogonal_subset;
    }

    PointCloud EstimateNormals(PointCloud pc, double knn_rad, int32_t max_nn){
        if (!pc.hasNormals()){
            pc.estimateNormalsKNNInRadius(max_nn, knn_rad, true);
        }
        
        return BuildNormalsAndLambdas(pc, knn_rad);
    }

    PointCloud BuildNormalsAndLambdas(PointCloud const & pc, double knn_rad){
        auto pc_tree = std::make_unique<KDTree>(pc.points);
        auto points = pc.points;
        auto main_normals = pc.normals;

        std::vector<decltype(main_normals.col(0))> normals;
        std::vector<decltype(points.col(0))> new_points;

        for (Eigen::Index i = 0; i < points.cols(); ++i){
            auto point = points.col(i);
            std::vector<Eigen::Index> indices = metrics_utils::GetRadiusSearchIndices(*pc_tree, point, knn_rad);
            if (indices.size() <= 3) continue;
            Eigen::MatrixX3d tmp = metrics_utils::TransformPointIdxToMatrix(points, indices);
            Eigen::MatrixX3d cov_matrix = metrics_utils::FindCov(tmp);
            Eigen::VectorXd eigenvalues = cov_matrix.eigenvalues().real();
            std::sort(eigenvalues.begin(), eigenvalues.end());
            if (100 * eigenvalues[0] < eigenvalues[1]){
                normals.push_back(main_normals.col(i));
                new_points.push_back(point);
            }
        }

        // TODO: Simplify transformations

        // new_points[0] shape (3, 1)  i.e. (x, y, z)^T
        cilantro::VectorSet3d pc_points(3, new_points.size());
        cilantro::VectorSet3d pc_normals(3, normals.size());
        
        for (Eigen::Index i = 0; i < new_points.size(); ++i){
            pc_points.col(i) = new_points[i];
            pc_normals.col(i) = normals[i];
        }
        
        auto new_pc = PointCloud();
    
        new_pc.points = pc_points;
        new_pc.normals = pc_normals;

        return new_pc;
    }
}
