#ifndef MAP_METRICS_CLUSTERING_H
#define MAP_METRICS_CLUSTERING_H

#include <vector>
#include <cilantro/utilities/point_cloud.hpp>
#include <ap.h>

namespace clustering{

    class ClusterLabels{
     public:
        explicit ClusterLabels(alglib::integer_1d_array const & labels, Eigen::Index cluster_number);
     private:
        alglib::integer_1d_array labels_;
        Eigen::Index cluster_number_;
    };

    ClusterLabels ClusterizeAHC(Eigen::Ref<const cilantro::VectorSet3d> const points, 
                                   double distance_treshold);

}

#endif // MAP_METRICS_CLUSTERING_H