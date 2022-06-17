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
#ifndef MAP_METRICS_METRICS_H
#define MAP_METRICS_METRICS_H

#include "metrics_utils.h"
#include "config.h"

namespace metrics{
    using PointCloud = cilantro::PointCloud3d;
    using KDTree = cilantro::KDTree3d<>;

    // Available Base metrics are:
    // MPV -- Mean Plane Variance
    // MME -- Mean Map Entropy
    double ComputeBaseMetrics(
            std::vector<cilantro::VectorSet3d> const & pcs_points,
            std::vector<Eigen::Matrix4d> const & ts,
            config::CustomConfig config = config::LidarConfig(),
            std::optional<double> (*algorithm)
            (cilantro::VectorSet3d const & points,
             std::vector<Eigen::Index> const & indices) = metrics_utils::metrics_algorithm::ComputeEigenvalues
            );

    double ComputeOrthogonalMetrics(
            std::vector<cilantro::VectorSet3d> const & pcs_points,
            std::vector<Eigen::Matrix4d> const & ts, 
            config::CustomConfig config = config::LidarConfig(),
            std::optional<double> (*algorithm)
            (cilantro::VectorSet3d const & points,
             std::vector<Eigen::Index> const & indices) = metrics_utils::metrics_algorithm::ComputeEigenvalues,
            std::vector<cilantro::VectorSet3d> const & orth_subset = std::vector<cilantro::VectorSet3d>()
            );

    // MPV. Mean Plane Variance
    double GetMPV(std::vector<cilantro::VectorSet3d> const & pcs_points, std::vector<Eigen::Matrix4d> const & ts,
               config::CustomConfig config = config::LidarConfig());

    // MME. Mean Map Entropy
    double GetMME(std::vector<cilantro::VectorSet3d> const & pcs_points, std::vector<Eigen::Matrix4d> const & ts,
               config::CustomConfig config = config::LidarConfig());

    // MOM. Mutually Orthogonal Metric 
    double GetMOM(std::vector<cilantro::VectorSet3d> const & pcs_points, std::vector<Eigen::Matrix4d> const & ts,
               config::CustomConfig config = config::LidarConfig(), std::vector<cilantro::VectorSet3d> const & orth_subset = {});
} // namespace metrics

#endif //MAP_METRICS_METRICS_H
