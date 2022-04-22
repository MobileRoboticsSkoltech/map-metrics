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

} // namespace config

#endif // MAP_METRICS_CONFIG_H
