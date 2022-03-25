#ifndef MAP_METRICS_CONFIG_H
#define MAP_METRICS_CONFIG_H

namespace config{

    class CustomConfig{
    public:
        CustomConfig(int min_knn = 5, double knn_rad = 1.0, 
                    int max_nn = 30, int min_clust_size = 5): 
                    min_knn(min_knn), knn_rad(knn_rad), 
                    max_nn(max_nn), min_clust_size(min_clust_size) {}

        const int min_knn;
        const double knn_rad; 
        const int max_nn;
        const int min_clust_size;

    };

    class LidarConfig: public CustomConfig{
    public:
        LidarConfig(): CustomConfig() {}
    };

} // namespace config

#endif // MAP_METRICS_CONFIG_H
