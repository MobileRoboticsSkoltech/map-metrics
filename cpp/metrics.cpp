//
// Created by achains on 07.12.2021.
//

#include "metrics_utils.h"
#include "metrics.h"

namespace metrics {
    double compute_base_metric(
            std::vector<cilantro::VectorSet3d> const &pcs_points,
            std::vector<Eigen::Matrix4d> const &ts,
            int min_knn,
            double knn_radius,
            std::optional<double> (*algorithm)
                    (cilantro::VectorSet3d const &points,
                     std::vector<unsigned long> const &indices)) {

        std::vector<cilantro::PointCloud3d> pcs(pcs_points.size());
        for (size_t i = 0; i < pcs_points.size(); ++i){
            pcs[i] = cilantro::PointCloud3d(pcs_points[i]);
        }

        metrics_utils::PointCloud pc_map = metrics_utils::aggregate_map(pcs, ts);

        metrics_utils::KDTree map_tree(pc_map.points);
        cilantro::VectorSet3d points = pc_map.points;

        std::vector<double> metric;
        for (Eigen::Index i = 0; i < points.cols(); ++i) {
            std::vector<unsigned long> indices = metrics_utils::get_radius_search_indices(map_tree,
                                                                           points.col(i), knn_radius);
            if (indices.size() > min_knn) {
                std::optional<double> result = algorithm(points, indices);
                if (result.has_value())
                    metric.push_back(result.value());
            }
        }

        return (metric.empty() ? 0 :
        std::reduce(metric.begin(), metric.end()) / static_cast<double>(metric.size()));
    }

    double mpv(std::vector<cilantro::VectorSet3d> const & pcs_points, std::vector<Eigen::Matrix4d> const & ts,
               int min_knn, double knn_rad){
        return compute_base_metric(pcs_points, ts, min_knn, knn_rad,
                                   &metrics_utils::metrics_algorithm::compute_eigenvalues);
    }

    double mme(std::vector<cilantro::VectorSet3d> const & pcs_points, std::vector<Eigen::Matrix4d> const & ts,
               int min_knn, double knn_rad){
        return compute_base_metric(pcs_points, ts, min_knn, knn_rad,
                                   &metrics_utils::metrics_algorithm::compute_entropy);
    }
} // namespace metrics