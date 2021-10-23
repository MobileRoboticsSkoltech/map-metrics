#ifndef MAP_METRICS_METRICS_H
#define MAP_METRICS_METRICS_H

#include <algorithm>
#include <numeric>
#include <vector>

#include "open3d/Open3D.h"
#include "Eigen/Core"
#include "Eigen/Eigenvalues"

namespace map_metrics {
    using open3d::geometry::PointCloud;
    using open3d::geometry::KDTreeFlann;

    struct SearchRadiusData{
        int k;
        std::vector<int> indices;
        std::vector<double> distance2;
    };

    Eigen::MatrixX3d cov(Eigen::MatrixX3d const & M);

    Eigen::MatrixX3d points_idx_to_matrix(std::vector<Eigen::Vector3d> const & points, std::vector<int> const & idx);

    PointCloud aggregate_map(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts);

    SearchRadiusData search_radius_vector_3d(KDTreeFlann const & tree, Eigen::Vector3d const & query, double radius);
    
}

#endif //MAP_METRICS_METRICS_H