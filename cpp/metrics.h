//
// Created by achains on 18.11.2021.
//

#ifndef MAP_METRICS_METRICS_H
#define MAP_METRICS_METRICS_H

#include <cilantro/core/kd_tree.hpp>
#include <cilantro/utilities/point_cloud.hpp>

namespace map_metrics {

    using PointCloud = cilantro::PointCloud3d;
    using KDTree = cilantro::KDTree3d<>;

    Eigen::MatrixX3d cov(Eigen::MatrixX3d const & M);

    // Points are represented like
    //    x1   ...   xn
    //  ( y1 ) ... ( yn )
    //    z1   ...   zn
    Eigen::MatrixX3d points_idx_to_matrix(cilantro::VectorSet3d const & points, std::vector<unsigned long> const & idx);

    PointCloud aggregate_map(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts);

    std::vector<unsigned long> get_radius_search_indices(KDTree const & tree,
                                                       Eigen::Vector3d const & query, double radius);

    double mpv(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts);

    double mme(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts);

}

#endif //MAP_METRICS_METRICS_H
