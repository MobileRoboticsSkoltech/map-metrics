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
#ifndef MAP_METRICS_CONFIG_H
#define MAP_METRICS_CONFIG_H

#include <cstdint>

namespace config{

    class CustomConfig{
    public:
        CustomConfig(int32_t min_knn = 5, double knn_rad = 1.0, 
                    int32_t max_nn = 30, int32_t min_clust_size = 5): 
                    min_knn(min_knn), knn_rad(knn_rad), 
                    max_nn(max_nn), min_clust_size(min_clust_size) {}

        const int32_t min_knn;
        const double knn_rad; 
        const int32_t max_nn;
        const int32_t min_clust_size;

    };

    class LidarConfig: public CustomConfig{
    public:
        LidarConfig(): CustomConfig() {}
    };

    class DepthConfig: public CustomConfig{
    public:
        DepthConfig(): CustomConfig(5, 0.2) {}
    };

} // namespace config

#endif // MAP_METRICS_CONFIG_H
