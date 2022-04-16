//
// Created by achains on 07.12.2021.
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
               config::CustomConfig config = config::LidarConfig());
} // namespace metrics

#endif //MAP_METRICS_METRICS_H
