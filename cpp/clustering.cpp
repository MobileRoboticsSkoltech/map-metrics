#include "clustering.h"

#include <dataanalysis.h>

namespace clustering{
    ClusterMeans::ClusterMeans(alglib::integer_1d_array const & labels, Eigen::Index cluster_number) 
        : cluster_number_(cluster_number), labels_(Eigen::VectorXd(labels.length())),
          cluster_means_(), cluster_idx_(){
            for (Eigen::Index i = 0; i < labels.length(); ++i) labels_[i] = labels[i];
    }

    void ClusterMeans::filterClusters(Eigen::Ref<const cilantro::VectorSet3d> const points, int min_clust_size){
        cluster_idx_.resize(cluster_number_);
        cluster_means_.resize(cluster_number_);

        for (Eigen::Index i = 0; i < cluster_number_; ++i){
            // ...
        }
    }
    
    
    ClusterMeans ClusterizeAHC(Eigen::Ref<const cilantro::VectorSet3d> const points,
                                   double distance_treshold){
        // Define clusterizer
        alglib::clusterizerstate s;
        alglib::ahcreport rep;
        alglib::clusterizercreate(s);

        // Fit points data
        alglib::real_2d_array xy;
        xy.setcontent(points.cols(), 3, points.data());
        constexpr int disttype = 2;

        alglib::integer_1d_array labels;
        alglib::integer_1d_array cz;
        alglib::ae_int_t number_of_clusters;

        // Clustering
        alglib::clusterizersetpoints(s, xy, disttype);
        alglib::clusterizerrunahc(s, rep);

        // Get top clusters from AHC tree which are separated by distance treshold
        alglib::clusterizerseparatedbydist(rep, distance_treshold, number_of_clusters, labels, cz);

        return ClusterMeans(labels, number_of_clusters);
    }

}