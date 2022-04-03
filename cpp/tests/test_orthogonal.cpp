#include <gtest/gtest.h>

#include "parse_utils.h"

#include "metrics.h"
#include "orth_extract.h"

#include <iostream>
#include <dataanalysis.h>


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

    // Define clusterizer
    alglib::clusterizerstate s;
    alglib::ahcreport rep;
    alglib::clusterizercreate(s);

    // Fit Data
    alglib::real_2d_array xy;
    xy.setcontent(normals.cols(), 3, normals.data());
    constexpr int disttype = 2;

    // Labels variable definition
    alglib::integer_1d_array labels;
    alglib::integer_1d_array cz;
    alglib::ae_int_t number_of_clusters;

    // Clustering
    alglib::clusterizersetpoints(s, xy, disttype);
    alglib::clusterizerrunahc(s, rep);
    
    // top clusters from  hierarchical  clusterization  tree which are separated by distance >R
    // Each  of  the clusters stand on its own, no heirarchy
    alglib::clusterizerseparatedbydist(rep, 1e-1, number_of_clusters, labels, cz);

    //std::cout << labels.tostring().c_str() << std::endl;
    
}
