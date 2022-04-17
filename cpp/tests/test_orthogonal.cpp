#include <gtest/gtest.h>

#include "parse_utils.h"

#include "max_clique.h"
#include "metrics.h"
#include "orth_extract.h"
#include "clustering.h"

#include <iostream>
#include <dataanalysis.h>

#ifdef TEST_ORTHOGONAL_METRICS

TEST(Orthogonal, FilterClusters){
    std::vector<cilantro::VectorSet3d> PCs = parse_utils::depth_parse::GetDepthPCs("data/depth/pcs");
    auto pc = cilantro::PointCloud3d(PCs[2]);  // PC-0091.ply

    // Note: cut_pc gives same result as Python version
    auto cut_pc = orth_extract::EstimateNormals(pc, 0.2, 30);
    auto normals = cut_pc.normals;

    clustering::ClusterMeans cluster_means = clustering::ClusterizeAHC(normals, 1e-1);
    cluster_means.filterClusters(normals, 5);

    auto cluster_means_pc = cilantro::PointCloud3d();
    cluster_means_pc.points = cluster_means.getMeans();

    cluster_means_pc.toPLYFile("cluster_means_pc_cpp_2.ply");
}

TEST(Orthogonal, CliqueEstimation){
    std::vector<cilantro::VectorSet3d> PCs = parse_utils::depth_parse::GetDepthPCs("data/depth/pcs");
    auto pc = cilantro::PointCloud3d(PCs[2]);  // PC-0091.ply

    // Note: cut_pc gives same result as Python version
    auto cut_pc = orth_extract::EstimateNormals(pc, 0.2, 30);
    auto normals = cut_pc.normals;

    clustering::ClusterMeans cluster_means = clustering::ClusterizeAHC(normals, 1e-1);
    cluster_means.filterClusters(normals, 5);

    max_clique::FindMaxClique(cluster_means, 1e-1);
}

TEST(Orthogonal, OrthSubsetEstimation){
   std::vector<cilantro::VectorSet3d> PCs = parse_utils::depth_parse::GetDepthPCs("data/depth/pcs");
   auto pc = cilantro::PointCloud3d(PCs[2]);  // PC-0091.ply

   auto orth_subset = orth_extract::ExtractOrthogonalSubset(pc, config::CustomConfig(5, 0.2, 30, 5));

   cilantro::PointCloud3d orth_pc;
   for (auto const & points : orth_subset){
       orth_pc.append(cilantro::PointCloud3d(points));
   }

   orth_pc.toPLYFile("orth_PC_cpp.ply");
}

TEST(Orthogonal, MOM){
    std::vector<cilantro::VectorSet3d> PCs = parse_utils::depth_parse::GetDepthPCs("data/depth/pcs");
    std::vector<Eigen::Matrix4d> poses = parse_utils::depth_parse::GetDepthPoses("data/depth/poses");
    std::vector<Eigen::Matrix4d> auto_poses = parse_utils::depth_parse::GetDepthPoses("data/depth/auto-poses");
    auto pc = cilantro::PointCloud3d(PCs[2]);  // PC-0091.ply
//    auto pivot_pose = poses[2];
//    cilantro::RigidTransform3d transform_mx(pivot_pose);

    auto orth_subset = orth_extract::ExtractOrthogonalSubset(pc, config::CustomConfig(5, 0.2, 30, 5));

//    // orth_subset rigid transformation
//    for (size_t i = 0; i < orth_subset.size(); ++i){
//        auto pc_transformed = cilantro::PointCloud3d(orth_subset[i]).transformed(transform_mx);
//        orth_subset[i] = pc_transformed.points;
//    }

    std::cout << "MOM-init: " <<
        metrics::ComputeOrthogonalMetrics(
                PCs,
                poses,
                config::CustomConfig(5, 0.2, 30, 5),
                &metrics_utils::metrics_algorithm::ComputeEigenvalues,
                orth_subset
                )
    << std::endl;

        std::cout << "MOM-auto: " <<
        metrics::ComputeOrthogonalMetrics(
                PCs,
                auto_poses,
                config::CustomConfig(5, 0.2, 30, 5),
                &metrics_utils::metrics_algorithm::ComputeEigenvalues,
                orth_subset
                )
    << std::endl;

}

#endif // TEST_ORTHOGONAL_METRICS