//
// Created by achains on 18.11.2021.
//

#ifndef MAP_METRICS_METRICS_H
#define MAP_METRICS_METRICS_H


#include <iostream>
#include <cilantro/core/kd_tree.hpp>
#include <cilantro/utilities/point_cloud.hpp>

namespace map_metrics {

    using PointCloud = cilantro::PointCloud3d;
    using KDTree = cilantro::KDTree3d<>;

    struct SearchRadiusData{
        int k;
        std::vector<int> indices;
        std::vector<double> distance2;
    };

    Eigen::MatrixX3d cov(Eigen::MatrixX3d const & M);

    Eigen::MatrixX3d points_idx_to_matrix(std::vector<Eigen::Vector3d> const & points, std::vector<int> const & idx);

    PointCloud aggregate_map(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts);

    SearchRadiusData search_radius_vector_3d(KDTree const & tree, Eigen::Vector3d const & query, double radius);

    double mpv(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts);

    double mme(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts);

}

#endif //MAP_METRICS_METRICS_H
