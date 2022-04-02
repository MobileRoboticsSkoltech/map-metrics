#include <gtest/gtest.h>

#include "parse_utils.h"

#include "metrics.h"
#include "orth_extract.h"

#include <iostream>


// TEST(Orthogonal, NormalEstimation){
//     std::vector<cilantro::VectorSet3d> PCs = GetPointClouds("data/kitti_00");
//     auto pc = cilantro::PointCloud3d(PCs[0]);

//     pc = orth_extract::EstimateNormals(pc, 1.0, 30);
    
//     ASSERT_TRUE(pc.hasNormals());

//     auto normals = pc.normals;
//     std::vector<double> n(normals.col(0).begin(), normals.col(0).end());
//     //std::cout << normals.col(0) << std::endl;
// }

TEST(Orthogonal, BuildNormalsAndLambdas){
    std::vector<cilantro::VectorSet3d> PCs = parse_utils::depth_parse::GetDepthPCs("data/depth/pcs");
    auto pc = cilantro::PointCloud3d(PCs[2]);  // PC-0091.ply

    auto cut_pc = orth_extract::EstimateNormals(pc, 0.2, 30);
    auto normals = cut_pc.normals;
    
    // std::cout << "Normals(0):\n" << cut_pc.normals.col(0) << '\n';
    // std::cout << "Normals(25003):\n" << cut_pc.normals.col(25003) << '\n';
}
