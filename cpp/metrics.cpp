//
// Created by achains on 18.11.2021.
//

#include "metrics.h"

#include <algorithm>
#include <numeric>

#include "Eigen/Core"
#include "Eigen/Dense"

namespace map_metrics{

    Eigen::MatrixX3d cov(Eigen::MatrixX3d const & M){
        Eigen::MatrixX3d centered = M.rowwise() - M.colwise().mean();
        return (centered.adjoint() * centered) / (static_cast<double>(M.rows()) - 1.0);
    }

    Eigen::MatrixX3d points_idx_to_matrix(cilantro::VectorSet3d const & points, std::vector<unsigned long> const & idx){
        Eigen::MatrixX3d matrix(idx.size(), 3);
        Eigen::Index row_idx = 0;
        for (const auto & i : idx){
            matrix.row(row_idx++) = points.col(static_cast<Eigen::Index>(i));
        }
        return matrix;
    }

    PointCloud aggregate_map(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts){
        assert (pcs.size() == ts.size());

        std::vector<Eigen::Matrix4d> inv_ts;
        Eigen::Matrix4d inv_elem = ts[0].inverse();
        for (Eigen::Matrix4d const & mx : ts){
            inv_ts.emplace_back(mx * inv_elem);
        }

        cilantro::VectorSet3d pc_map_points(3, 0);
        for (size_t i = 0; i < pcs.size(); ++i){
            cilantro::RigidTransform3d transform_mx(inv_ts[i]);
            cilantro::VectorSet3d transformed_points = pcs[i].transformed(transform_mx).points;

            Eigen::Index old_size = pc_map_points.cols();
            pc_map_points.conservativeResize(3, old_size + transformed_points.cols());
            for (Eigen::Index col_idx = 0; col_idx < transformed_points.cols(); ++col_idx){
                pc_map_points.col(old_size + col_idx) = transformed_points.col(col_idx);
            }
        }

        return PointCloud{pc_map_points};
    }

    std::vector<unsigned long> get_radius_search_indices(KDTree const & tree,
                                                       Eigen::Vector3d const & query, double radius){
        cilantro::NeighborSet<double> nn = tree.radiusSearch(query, radius);
        std::vector<unsigned long> indices(nn.size());
        for (size_t i = 0; i < nn.size(); ++i){
            indices[i] = nn[i].index;
        }

        return indices;
    }

    double mpv(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts){
        const int min_knn = 5;
        const double knn_rad = 1.0;

        PointCloud pc_map = aggregate_map(pcs, ts);

        KDTree map_tree(pc_map.points);
        cilantro::VectorSet3d points = pc_map.points;

        std::vector<double> metric;
        for (Eigen::Index i = 0; i < points.cols(); ++i){
            std::vector<unsigned long> indices = get_radius_search_indices(map_tree,
                                                                        points.col(i), knn_rad);
            if (indices.size() > min_knn){
                Eigen::MatrixX3d tmp = points_idx_to_matrix(points, indices);
                Eigen::MatrixX3d cov_matrix = cov(tmp);
                Eigen::VectorXd eigenvalues = cov_matrix.eigenvalues().real();
                metric.push_back(eigenvalues.minCoeff());
            }
        }
        return (metric.empty() ? 0 : std::reduce(metric.begin(), metric.end()) / static_cast<double>(metric.size()));
    }

    double mme(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts){
        const int min_knn = 5;
        const double knn_rad = 1.0;

        PointCloud pc_map = aggregate_map(pcs, ts);

        KDTree map_tree(pc_map.points);
        cilantro::VectorSet3d points = pc_map.points;

        std::vector<double> metric;
        for (Eigen::Index i = 0; i < points.cols(); ++i){
            std::vector<unsigned long> indices = get_radius_search_indices(map_tree,
                                                                           points.col(i), knn_rad);
            if (indices.size() > min_knn){
                Eigen::MatrixXd tmp = points_idx_to_matrix(points, indices);
                Eigen::MatrixXd cov_matrix = cov(tmp);
                double det =  (2. * M_PI * M_E * cov_matrix).determinant();
                if (det > 0)
                    metric.push_back(0.5 * std::log(det));
            }
        }
        return (metric.empty() ? 0 : std::reduce(metric.begin(), metric.end()) / static_cast<double>(metric.size()));
    }

}