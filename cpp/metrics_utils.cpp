//
// Created by achains on 18.11.2021.
//

#include "metrics_utils.h"

#include <algorithm>
#include <cmath>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/StdVector"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.718281828459045235360
#endif

namespace metrics_utils{

    Eigen::MatrixX3d FindCov(Eigen::MatrixX3d const & M){
        Eigen::MatrixX3d centered = M.rowwise() - M.colwise().mean();
        return (centered.adjoint() * centered) / (static_cast<double>(M.rows()) - 1.0);
    }

    Eigen::MatrixX3d TransformPointIdxToMatrix(cilantro::VectorSet3d const & points, std::vector<int> const & idx){
        Eigen::MatrixX3d matrix(idx.size(), 3);
        Eigen::Index row_idx = 0;
        for (const auto & i : idx){
            matrix.row(row_idx++) = points.col(static_cast<Eigen::Index>(i));
        }
        return matrix;
    }

    PointCloud AggregateMap(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts){
        assert (pcs.size() == ts.size());

        std::vector<Eigen::Matrix4d> inv_ts;
        Eigen::Matrix4d inv_elem = ts[0].inverse();
        for (const Eigen::Ref<const Eigen::Matrix4d> mx : ts){
            inv_ts.push_back(inv_elem * mx);
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

    std::vector<int> GetRadiusSearchIndices(KDTree const & tree,
                                                       Eigen::Vector3d const & query, double radius){
        cilantro::NeighborSet<double> nn = tree.radiusSearch(query, radius);
        std::vector<int> indices(nn.size());
        for (size_t i = 0; i < nn.size(); ++i){
            indices[i] = static_cast<int>(nn[i].index);
        }

        return indices;
    }

    namespace metrics_algorithm{
        std::optional<double> ComputeEigenvalues(cilantro::VectorSet3d const & points,
                                   std::vector<int> const & indices){
            Eigen::MatrixX3d tmp = TransformPointIdxToMatrix(points, indices);
            Eigen::MatrixX3d cov_matrix = FindCov(tmp);
            Eigen::VectorXd eigenvalues = cov_matrix.eigenvalues().real();
            return eigenvalues.minCoeff();
        }

        std::optional<double> ComputeEntropy(cilantro::VectorSet3d const & points,
                                              std::vector<int> const & indices){
            Eigen::MatrixXd tmp = TransformPointIdxToMatrix(points, indices);
            Eigen::MatrixXd cov_matrix = FindCov(tmp);
            double det =  (2. * M_PI * M_E * cov_matrix).determinant();
            if (det > 0)
                return 0.5 * std::log(det);
            return {};
        }
    } // namespace metrics_algorithm
} // namespace metrics_utils