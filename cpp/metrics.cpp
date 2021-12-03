//
// Created by achains on 18.11.2021.
//

#include "metrics.h"

#include <algorithm>
#include <numeric>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Geometry"

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

        PointCloud pc_map = PointCloud();
        for (size_t i = 0; i < pcs.size(); ++i){
            // TODO: Fix transform
            Eigen::Transform <double, 3, Eigen::Affine> t(inv_ts[i]);
            pc_map.points += pcs[i].transformed(t).points;
        }

        return pc_map;
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

//    double mme(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts){
//        const int min_knn = 5;
//        const double knn_rad = 1.0;
//
//        PointCloud pc_map = aggregate_map(pcs, ts);
//
//        KDTree map_tree(pc_map.points);
//
//        cilantro::VectorSet3d points = pc_map.points;
//
//        std::vector<double> metric;
//        for (const auto& point: points){
//            auto tpl = search_radius_vector_3d(map_tree, point, knn_rad);
//            if (tpl.indices.size() > min_knn){
//                Eigen::MatrixXd tmp = points_idx_to_matrix(points, tpl.indices);
//                Eigen::MatrixXd cov_matrix = cov(tmp);
//                double det =  (2. * M_PI * M_E * cov_matrix).determinant();
//                if (det > 0)
//                    metric.push_back(0.5 * std::log(det));
//            }
//        }
//        return (metric.empty() ? 0 : std::reduce(metric.begin(), metric.end()) / static_cast<double>(metric.size()));
//    }

}