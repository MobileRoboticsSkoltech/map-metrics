#include <gtest/gtest.h>

#include "cilantro/3rd_party/tinyply/tinyply.h"

#include "parse_utils.h"
#include "pc_utils.h"

#include "metrics.h"
#include <iostream>

namespace depth_parse{
    std::vector<cilantro::VectorSet3d> GetDepthPCs(char const * path){
        std::vector<std::filesystem::path> PC_file_names;
        for (const auto & entry : std::filesystem::directory_iterator(path)){
            PC_file_names.push_back(entry.path());
        }
        std::sort(PC_file_names.begin(), PC_file_names.end());

        std::vector<cilantro::VectorSet3d> pcs;
        for (auto const & filename : PC_file_names){
            cilantro::PointCloud3d pc;
            pc.fromPLYFile(filename);
            pcs.emplace_back(pc.points);
        }

        return pcs;
    }

    std::vector<cilantro::PointCloud3d> GetDepthPCsClouds(char const * path){
        std::vector<std::filesystem::path> PC_file_names;
        for (const auto & entry : std::filesystem::directory_iterator(path)){
            PC_file_names.push_back(entry.path());
        }
        std::sort(PC_file_names.begin(), PC_file_names.end());

        std::vector<cilantro::PointCloud3d> pcs;
        for (auto const & filename : PC_file_names){
            cilantro::PointCloud3d pc;
            pc.fromPLYFile(filename);
            pcs.emplace_back(pc);
        }

        return pcs;
    }

    std::vector<Eigen::Matrix4d> GetDepthPoses(char const * path){
        std::vector<std::filesystem::path> PC_file_names;
        for (const auto & entry : std::filesystem::directory_iterator(path)){
            PC_file_names.push_back(entry.path());
        }
        std::sort(PC_file_names.begin(), PC_file_names.end());

        std::vector<Eigen::Matrix4d> poses;
        for (auto const & filename : PC_file_names){
            poses.emplace_back(
                pc_utils::VectorToEigenMatrix(parse_utils::ParseTrajectory(filename))
            );
        }

        return poses;
    }
}

TEST(DepthTest, MME){
    std::vector<cilantro::VectorSet3d> pcs = depth_parse::GetDepthPCs("data/depth/pcs");
    //auto pcs = depth_parse::GetDepthPCsClouds("data/depth/pcs");
    std::vector<Eigen::Matrix4d> poses = depth_parse::GetDepthPoses("data/depth/poses");

    ASSERT_EQ(pcs.size(), 28);
    ASSERT_EQ(pcs.size(), poses.size());

    // Poses equal to Python result
    // PCs equal to Python result

    //auto pc_map = metrics_utils::AggregateMap(pcs, poses);

    std::cout << "MME: " << metrics::GetMME(pcs, poses, config::CustomConfig(5, 0.2, 30, 5)) << std::endl;


    // std::cout << "MPV: " << metrics::GetMPV(pcs, poses, config::CustomConfig(5, 0.2, 30, 5)) << std::endl;
}
