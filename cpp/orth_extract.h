// Copyright (c) 2022, Skolkovo Institute of Science and Technology (Skoltech)
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//
//  Created on: May 20, 2022
//       Author: Arthur Saliou
//               arthur.salio@gmail.com
//
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
