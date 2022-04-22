#ifndef MAP_METRICS_CLUSTERING_H
#define MAP_METRICS_CLUSTERING_H

#include <cstdint>
#include <vector>
#include <cilantro/core/data_containers.hpp>
#include <alglib/ap.h>

namespace clustering{

    class ClusterMeans{
     public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        explicit ClusterMeans(alglib::integer_1d_array const & labels, Eigen::Index cluster_number);

        void filterClusters(Eigen::Ref<const cilantro::VectorSet3d> const points, int32_t min_clust_size);
        
        const cilantro::VectorSet3d & getMeans() const;

        const Eigen::VectorXi & getIdx() const;

        const Eigen::VectorXi & getLabels() const;

        Eigen::Index getClusterNumber() const;

     private:
        Eigen::VectorXi labels_;
        cilantro::VectorSet3d cluster_means_;
        Eigen::VectorXi cluster_idx_;
        Eigen::Index cluster_number_;
    };

    ClusterMeans ClusterizeAHC(Eigen::Ref<const cilantro::VectorSet3d> const points, 
                                   double distance_treshold);

}

#endif // MAP_METRICS_CLUSTERING_H