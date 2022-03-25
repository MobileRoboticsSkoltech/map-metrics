//
// Created by achains on 07.12.2021.
//

#include "metrics_utils.h"
#include "metrics.h"

#include <numeric>

namespace metrics {
    double ComputeBaseMetrics(
            std::vector<cilantro::VectorSet3d> const &pcs_points,
            std::vector<Eigen::Matrix4d> const &ts,
            config::CustomConfig config,
            std::optional<double> (*algorithm)
                    (cilantro::VectorSet3d const &points,
                     std::vector<int> const &indices)) {

        std::vector<cilantro::PointCloud3d> pcs(pcs_points.size());
        for (size_t i = 0; i < pcs_points.size(); ++i){
            pcs[i] = cilantro::PointCloud3d(pcs_points[i]);
        }

        metrics_utils::PointCloud pc_map = metrics_utils::AggregateMap(pcs, ts);

        metrics_utils::KDTree map_tree(pc_map.points);
        cilantro::VectorSet3d points = pc_map.points;

        std::vector<double> metric;
        for (Eigen::Index i = 0; i < points.cols(); ++i) {
            std::vector<int> indices = metrics_utils::GetRadiusSearchIndices(map_tree,
                                                                           points.col(i), config.knn_rad);
            if (indices.size() > config.min_knn) {
                std::optional<double> result = algorithm(points, indices);
                if (result.has_value())
                    metric.push_back(result.value());
            }
        }

        return (metric.empty() ? 0 :
        std::accumulate(metric.begin(), metric.end(), 0.0) / static_cast<double>(metric.size()));
    }

    double ComputeOrthogonalMetrics(
            std::vector<cilantro::VectorSet3d> const & pcs_points,
            std::vector<Eigen::Matrix4d> const & ts, 
            config::CustomConfig config,
            std::optional<double> (*algorithm)
            (cilantro::VectorSet3d const & points,
             std::vector<int> const & indices)){
        
        std::vector<cilantro::PointCloud3d> pcs(pcs_points.size());
        for (size_t i = 0; i < pcs_points.size(); ++i){
            pcs[i] = cilantro::PointCloud3d(pcs_points[i]);
        }

        metrics_utils::PointCloud pc_map = metrics_utils::AggregateMap(pcs, ts);

        metrics_utils::KDTree map_tree(pc_map.points);
        cilantro::VectorSet3d points = pc_map.points;

        // auto orth_list = ExtractOrthogonalSubset
    }

    double GetMPV(std::vector<cilantro::VectorSet3d> const & pcs_points, std::vector<Eigen::Matrix4d> const & ts,
               config::CustomConfig config){
        return ComputeBaseMetrics(pcs_points, ts, config,
                                   &metrics_utils::metrics_algorithm::ComputeEigenvalues);
    }

    double GetMME(std::vector<cilantro::VectorSet3d> const & pcs_points, std::vector<Eigen::Matrix4d> const & ts,
               config::CustomConfig config){
        return ComputeBaseMetrics(pcs_points, ts, config,
                                   &metrics_utils::metrics_algorithm::ComputeEntropy);
    }

    
} // namespace metrics
