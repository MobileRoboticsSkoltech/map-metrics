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

    Eigen::MatrixX3d points_idx_to_matrix(std::vector<Eigen::Vector3d> const & points, std::vector<int> const & idx){
        Eigen::MatrixX3d matrix(idx.size(), 3);
        Eigen::Index row_idx = 0;
        for (const auto & i : idx){
            matrix.row(row_idx++) = points[i];
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
            PointCloud tmp = pcs[i];
            pc_map += tmp.Transform(inv_ts[i]);
        }

        return pc_map;
    }

    // https://github.com/isl-org/Open3D/blob/master/cpp/pybind/geometry/kdtreeflann.cpp#L176
    SearchRadiusData search_radius_vector_3d(KDTreeFlann const & tree, Eigen::Vector3d const & query, double radius){
        std::vector<int> indices;
        std::vector<double> distance2;
        int k = tree.SearchRadius(query, radius, indices, distance2);

        if (k < 0)
            throw std::runtime_error(
                    "search_radius_vector_3d() error!");

        return SearchRadiusData({k, indices, distance2});
    }

    double mpv(std::vector<PointCloud> const & pcs, std::vector<Eigen::Matrix4d> const & ts){
        const int min_knn = 5;
        const double knn_rad = 1.0;

        PointCloud pc_map = aggregate_map(pcs, ts);


        KDTreeFlann map_tree;
        map_tree.SetGeometry(pc_map);
        std::vector<Eigen::Vector3d> points = pc_map.points_;

        std::vector<double> metric;
        for (const auto& point: points){
            auto tpl = search_radius_vector_3d(map_tree, point, knn_rad);
            if (tpl.indices.size() > min_knn){
                Eigen::MatrixX3d tmp = points_idx_to_matrix(points, tpl.indices);
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

        KDTreeFlann map_tree;
        map_tree.SetGeometry(pc_map);
        std::vector<Eigen::Vector3d> points = pc_map.points_;

        std::vector<double> metric;
        for (const auto& point: points){
            auto tpl = search_radius_vector_3d(map_tree, point, knn_rad);
            if (tpl.indices.size() > min_knn){
                Eigen::MatrixXd tmp = points_idx_to_matrix(points, tpl.indices);
                Eigen::MatrixXd cov_matrix = cov(tmp);
                double det =  (2. * M_PI * M_E * cov_matrix).determinant();
                if (det > 0)
                    metric.push_back(0.5 * std::log(det));
            }
        }
        return (metric.empty() ? 0 : std::reduce(metric.begin(), metric.end()) / static_cast<double>(metric.size()));
    }

}