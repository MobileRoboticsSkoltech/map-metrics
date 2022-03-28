// #include <gtest/gtest.h>

// #include "pc_utils.h"
// #include "parse_utils.h"

// #include "metrics.h"
// #include "orth_extract.h"

// #include <iostream>

// std::vector<std::vector<double>> GetCalibrationMatrix(){
//     return
//             {
//                     {4.276802385584e-04, -9.999672484946e-01, -8.084491683471e-03, -1.198459927713e-02},
//                     {-7.210626507497e-03, 8.081198471645e-03, -9.999413164504e-01, -5.403984729748e-02},
//                     {9.999738645903e-01, 4.859485810390e-04, -7.206933692422e-03, -2.921968648686e-01},
//                     {0, 0, 0, 1}
//             };
// }

// std::vector<Eigen::Matrix4d> GetCalibratedTrajectory(std::filesystem::path const & path){
//     std::vector<std::vector<double>> trajectory = parse_utils::ParseTrajectory("data/00.txt");

//     return pc_utils::CalibrateTrajectory(GetCalibrationMatrix(), trajectory);
// }

// std::vector<cilantro::VectorSet3d> GetPointClouds(std::filesystem::path const & path){
//     std::vector<std::vector<double>> PCs_vector = parse_utils::GetPointClouds(path);
//     std::vector<cilantro::VectorSet3d> PCs_points;

//     for (auto & PC_vector : PCs_vector)
//     {
//         cilantro::VectorSet3d points = pc_utils::VectorToPointCloudPoints(PC_vector);
//         PCs_points.push_back(points);
//     }

//     return PCs_points;
// }

// // TEST(Orthogonal, NormalEstimation){
// //     std::vector<cilantro::VectorSet3d> PCs = GetPointClouds("data/kitti_00");
// //     auto pc = cilantro::PointCloud3d(PCs[0]);

// //     pc = orth_extract::EstimateNormals(pc, 1.0, 30);
    
// //     ASSERT_TRUE(pc.hasNormals());

// //     auto normals = pc.normals;
// //     std::vector<double> n(normals.col(0).begin(), normals.col(0).end());
// //     //std::cout << normals.col(0) << std::endl;
// // }

// // TEST(Orthogonal, BuildNormalsAndLambdas){
// //     std::vector<cilantro::VectorSet3d> PCs = GetPointClouds("data/kitti_00");
// //     auto pc = cilantro::PointCloud3d(PCs[0]);

// //     auto cut_pc = orth_extract::EstimateNormals(pc, 1.0, 30);

// //     std::cout << "Cut PC points/normals size: " << cut_pc.points.cols() << '\n';
// //     std::cout << "Normals(0):\n" << cut_pc.normals.col(0) << '\n';
// //     std::cout << "Normals(25003):\n" << cut_pc.normals.col(25003) << '\n';
// // }
