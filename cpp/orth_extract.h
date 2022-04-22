#ifndef MAP_METRICS_ORTH_EXTRACT_H
#define MAP_METRICS_ORTH_EXTRACT_H

#include "config.h"
#include "clustering.h"

#include <cstdint>
#include <cilantro/utilities/point_cloud.hpp>

namespace orth_extract{
    using PointCloud = cilantro::PointCloud3d;

    std::vector<cilantro::VectorSet3d> ExtractOrthogonalSubset(const PointCloud& pivot_pc,
                                    config::CustomConfig config = config::LidarConfig(),
                                    double eps = 1e-1);

    PointCloud EstimateNormals(PointCloud pc, double knn_rad, int32_t max_nn);

    PointCloud BuildNormalsAndLambdas(PointCloud const & pc, double knn_rad);
}

#endif //MAP_METRICS_ORTH_EXTRACT_H
