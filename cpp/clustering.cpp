#include "clustering.h"

#include <alglib/dataanalysis.h>

namespace clustering{
    ClusterMeans::ClusterMeans(alglib::integer_1d_array const & labels, Eigen::Index cluster_number) 
        : cluster_number_(cluster_number), labels_(Eigen::VectorXi(labels.length())),
          cluster_means_(), cluster_idx_(){
            for (Eigen::Index i = 0; i < labels.length(); ++i) labels_[i] = labels[i];
    }

    const cilantro::VectorSet3d & ClusterMeans::getMeans() const{
        return this->cluster_means_;
    }

    const Eigen::VectorXi & ClusterMeans::getIdx() const{
        return this->cluster_idx_;
    }

    const Eigen::VectorXi & ClusterMeans::getLabels() const{
        return this->labels_;
    }

    Eigen::Index ClusterMeans::getClusterNumber() const{
        return this->cluster_number_;
    }

    void ClusterMeans::filterClusters(Eigen::Ref<const cilantro::VectorSet3d> const points, int32_t min_clust_size){
        cluster_idx_.resize(cluster_number_);
        cluster_means_.resize(3, cluster_number_);

        Eigen::Index big_cluster_size = 0;

        for (Eigen::Index i = 0; i < cluster_number_; ++i){
            auto idx = (labels_.array() == i);
            if (idx.count() > min_clust_size){
                // TODO: Clean up! Find better solution for 'np.where'
                cilantro::VectorSet3d normals_at_idx(3, idx.count());
                Eigen::Index normals_at_idx_size = 0;
                for (Eigen::Index col_idx = 0; col_idx < idx.size(); ++col_idx){
                    if (idx[col_idx]) normals_at_idx.col(normals_at_idx_size++) = points.col(col_idx);
                }

                cluster_means_.col(big_cluster_size) = normals_at_idx.rowwise().mean().normalized();
                cluster_idx_[big_cluster_size] = i;
                ++big_cluster_size;
            }
        }

        // Shrink to fit
        cluster_idx_.conservativeResize(big_cluster_size);
        cluster_means_.conservativeResize(3, big_cluster_size);
        this->cluster_number_ = big_cluster_size;
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
        constexpr int32_t disttype = 2;

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