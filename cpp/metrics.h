//
// Created by achains on 07.12.2021.
//

#ifndef MAP_METRICS_METRICS_H
#define MAP_METRICS_METRICS_H

#include "metrics_utils.h"

namespace metrics{

    // Available Base metrics are:
    // MPV -- Mean Plane Variance
    // MME -- Mean Map Entropy
    double compute_base_metric(
            std::vector<cilantro::VectorSet3d> const & pcs_points,
            std::vector<Eigen::Matrix4d> const & ts,
            int min_knn = 5,
            double knn_radius = 1.0,
            std::optional<double> (*algorithm)
            (cilantro::VectorSet3d const & points,
             std::vector<int> const & indices) = metrics_utils::metrics_algorithm::compute_eigenvalues
            );

    // MPV. Mean Plane Variance
    double mpv(std::vector<cilantro::VectorSet3d> const & pcs_points, std::vector<Eigen::Matrix4d> const & ts,
               int min_knn = 5, double knn_rad = 1.0);

    // MME. Mean Map Entropy
    double mme(std::vector<cilantro::VectorSet3d> const & pcs_points, std::vector<Eigen::Matrix4d> const & ts,
               int min_knn = 5, double knn_rad = 1.0);
} // namespace metrics

#endif //MAP_METRICS_METRICS_H
