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

    Eigen::MatrixX3d TransformPointIdxToMatrix(cilantro::VectorSet3d const & points, std::vector<Eigen::Index> const & idx){
        Eigen::MatrixX3d matrix(idx.size(), 3);
        Eigen::Index row_idx = 0;
        for (const auto & i : idx){
            matrix.row(row_idx++) = points.col(i);
        }
        return matrix;
    }

    std::shared_ptr<PointCloud> AggregateMap(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts){
        assert (pcs.size() == ts.size());

        std::vector<Eigen::Matrix4d> inv_ts(ts.size());
        Eigen::Matrix4d inv_elem = ts[0].inverse();
        for (Eigen::Index i = 0; i < ts.size(); ++i){
            inv_ts[i] = inv_elem * ts[i];
        }

        auto pc_map = std::make_shared<PointCloud>();

        for (size_t i = 0; i < pcs.size(); ++i){
            cilantro::RigidTransform3d transform_mx(inv_ts[i]);
            PointCloud pc_transformed = pcs[i].transformed(transform_mx);
            pc_map->append(pc_transformed);
        }

        return pc_map;
    }

    std::vector<Eigen::Index> GetRadiusSearchIndices(KDTree const & tree,
                                                     const Eigen::Ref<const Eigen::Vector3d> query, double radius){
        cilantro::NeighborSet<double> nn = tree.radiusSearch(query, pow(radius, 2.0));
        std::vector<Eigen::Index> indices(nn.size());
        for (size_t i = 0; i < nn.size(); ++i){
            indices[i] = nn[i].index;
        }

        return indices;
    }

    namespace metrics_algorithm{
        std::optional<double> ComputeEigenvalues(cilantro::VectorSet3d const & points,
                                   std::vector<Eigen::Index> const & indices){
            Eigen::MatrixX3d tmp = TransformPointIdxToMatrix(points, indices);
            Eigen::MatrixX3d cov_matrix = FindCov(tmp);
            Eigen::VectorXd eigenvalues = cov_matrix.eigenvalues().real();
            return eigenvalues.minCoeff();
        }

        std::optional<double> ComputeEntropy(cilantro::VectorSet3d const & points,
                                              std::vector<Eigen::Index> const & indices){
            Eigen::MatrixXd tmp = TransformPointIdxToMatrix(points, indices);
            Eigen::MatrixXd cov_matrix = FindCov(tmp);
            double det =  (2. * M_PI * M_E * cov_matrix).determinant();
            if (det > 0)
                return 0.5 * std::log(det);
            return {};
        }
    } // namespace metrics_algorithm
} // namespace metrics_utils