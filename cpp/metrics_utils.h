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
#ifndef MAP_METRICS_METRICS_UTILS_H
#define MAP_METRICS_METRICS_UTILS_H

#include <cilantro/core/kd_tree.hpp>
#include <cilantro/utilities/point_cloud.hpp>

#include <optional>

namespace metrics_utils {

    using PointCloud = cilantro::PointCloud3d;
    using KDTree = cilantro::KDTree3d<>;

    Eigen::MatrixX3d FindCov(Eigen::MatrixX3d const & M);

    // Points are represented like
    //    x1   ...   xn
    //  ( y1 ) ... ( yn )
    //    z1   ...   zn
    Eigen::MatrixX3d TransformPointIdxToMatrix(cilantro::VectorSet3d const & points, std::vector<Eigen::Index> const & idx);

    std::shared_ptr<PointCloud> AggregateMap(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts);

    std::vector<Eigen::Index> GetRadiusSearchIndices(KDTree const & tree,
                                                       const Eigen::Ref<const Eigen::Vector3d> query, double radius);

    namespace metrics_algorithm{
        std::optional<double> ComputeEigenvalues(cilantro::VectorSet3d const & points,
                                                  std::vector<Eigen::Index> const & indices);

        std::optional<double> ComputeEntropy(cilantro::VectorSet3d const & points,
                                              std::vector<Eigen::Index> const & indices);
    } // namespace metrics_algorithm
} // namespace metrics_utils

#endif //MAP_METRICS_METRICS_UTILS_H
