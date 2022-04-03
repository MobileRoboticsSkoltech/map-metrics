#include "orth_extract.h"
#include "metrics_utils.h"
#include "clustering.h"

#include <cilantro/core/kd_tree.hpp>
#include <dataanalysis.h>

#include <iostream>

namespace orth_extract{
    using PointCloud = cilantro::PointCloud3d;
    using KDTree = cilantro::KDTree3d<>;

    void ExtractOrthogonalSubset(PointCloud pivot_pc, config::CustomConfig config, double eps){
        PointCloud cut_pc = BuildNormalsAndLambdas(pivot_pc, config.knn_rad);
        auto normals = cut_pc.normals;

        clustering::ClusterMeans labels = clustering::ClusterizeAHC(normals, 1e-1);
        labels.filterClusters(normals, config.min_clust_size);

        auto max_clique = FindMaxClique(labels, 1e-1);
    }

    PointCloud EstimateNormals(PointCloud pc, double knn_rad, int max_nn){
        if (!pc.hasNormals()){
            pc.estimateNormalsKNNInRadius(max_nn, knn_rad, true);
        }
        
        return BuildNormalsAndLambdas(pc, knn_rad);
    }

    PointCloud BuildNormalsAndLambdas(PointCloud const & pc, double knn_rad){
        auto pc_tree = std::make_unique<KDTree>(pc.points);
        auto points = pc.points;
        auto main_normals = pc.normals;

        std::vector<decltype(main_normals.col(0))> normals;
        std::vector<decltype(points.col(0))> new_points;

        for (Eigen::Index i = 0; i < points.cols(); ++i){
            auto point = points.col(i);
            std::vector<Eigen::Index> indices = metrics_utils::GetRadiusSearchIndices(*pc_tree, point, knn_rad);
            if (indices.size() <= 3) continue;
            Eigen::MatrixX3d tmp = metrics_utils::TransformPointIdxToMatrix(points, indices);
            Eigen::MatrixX3d cov_matrix = metrics_utils::FindCov(tmp);
            Eigen::VectorXd eigenvalues = cov_matrix.eigenvalues().real();
            std::sort(eigenvalues.begin(), eigenvalues.end());
            if (100 * eigenvalues[0] < eigenvalues[1]){
                normals.push_back(main_normals.col(i));
                new_points.push_back(point);
            }
        }

        // TODO: Simplify transformations

        // new_points[0] shape (3, 1)  i.e. (x, y, z)^T
        cilantro::VectorSet3d pc_points(3, new_points.size());
        cilantro::VectorSet3d pc_normals(3, normals.size());
        
        for (Eigen::Index i = 0; i < new_points.size(); ++i){
            pc_points.col(i) = new_points[i];
            pc_normals.col(i) = normals[i];
        }
        
        auto new_pc = PointCloud();
    
        new_pc.points = pc_points;
        new_pc.normals = pc_normals;

        return new_pc;
    }

    int FindMaxClique(clustering::ClusterMeans const & clusterizer, double eps){
        return 0;
    }
}