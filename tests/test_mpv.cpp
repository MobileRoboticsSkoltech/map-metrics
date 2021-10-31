//
// Created by achains on 30.10.2021.
//

#include <gtest/gtest.h>

#include "pc_utils.h"
#include "parse_utils.h"

#include "metrics.h"


TEST(Trajectory, TrajectoryCalibration){
    auto trajectory = parse_utils::ParseTrajectory("data/00.txt");
    std::vector<std::vector<double>> calibration_matrix =
            {
                    {4.276802385584e-04, -9.999672484946e-01, -8.084491683471e-03, -1.198459927713e-02},
                    {-7.210626507497e-03, 8.081198471645e-03, -9.999413164504e-01, -5.403984729748e-02},
                    {9.999738645903e-01, 4.859485810390e-04, -7.206933692422e-03, -2.921968648686e-01},
                    {0, 0, 0, 1}
            };

    std::vector<Eigen::Matrix4d> ts_gt = pc_utils::CalibrateTrajectory(calibration_matrix, trajectory);
    ASSERT_EQ(ts_gt.size(), 20);

    auto PCs_vector = parse_utils::GetPointClouds("data/kitti_00");
    std::vector<open3d::geometry::PointCloud> PCs;
    for (auto & PC_vector : PCs_vector)
    {
        auto PC = open3d::geometry::PointCloud();
        PC.points_ = pc_utils::VectorToPointCloudPoints(PC_vector);
        PCs.push_back(PC);
    }

    // PointCloud Data check
    // parse_utils::WritePCToFile("cpp_data.txt", PCs[0]);

    // Tj Data check
    parse_utils::WriteTjToFile("cpp_TJ_data.txt", ts_gt);
}
