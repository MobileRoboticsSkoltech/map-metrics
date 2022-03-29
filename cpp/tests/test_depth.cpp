#include <gtest/gtest.h>

#include "cilantro/3rd_party/tinyply/tinyply.h"

#include "parse_utils.h"
#include "pc_utils.h"

#include "metrics.h"
#include <iostream>

namespace depth_parse{
    std::vector<cilantro::VectorSet3d> GetDepthPCs(const char * path){
        std::vector<std::filesystem::path> pc_names;
        for (const auto & entry : std::filesystem::directory_iterator(path)){
            pc_names.push_back(entry.path());
        }
        std::sort(pc_names.begin(), pc_names.end());

        std::vector<cilantro::VectorSet3d> pcs;
        for (const auto & entry : pc_names){
            cilantro::PointCloud3d pc;
            pc.fromPLYFile(entry.string());
            pcs.push_back(pc.points);
        }

        return pcs;
    }

    std::vector<double> ParsePoses(std::filesystem::path path){
        std::ifstream input_file(path);
        std::string line;

        std::vector<double> poses;

        while(std::getline(input_file, line)){
            std::istringstream iss(line);
            std::string token;
            while(std::getline(iss, token, ' ')){
                poses.push_back(std::stod(token));
            }
        }

        return poses;
    }

    std::vector<Eigen::Matrix4d> GetDepthPoses(const char * path){
        std::vector<std::filesystem::path> poses_names;
        for (const auto & entry : std::filesystem::directory_iterator(path)){
            poses_names.push_back(entry.path());
        }
        std::sort(poses_names.begin(), poses_names.end());

        std::vector<Eigen::Matrix4d> poses; 
        for (const auto & entry : poses_names){
            auto p = ParsePoses(entry);
            Eigen::Matrix4d mx(p.data());
            poses.push_back(mx.transpose());
        }

        return poses;
    }
}

TEST(DepthTest, MME){
    std::vector<cilantro::VectorSet3d> pcs = depth_parse::GetDepthPCs("data/depth/pcs");
    // //auto pcs = depth_parse::GetDepthPCsClouds("data/depth/pcs");
    std::vector<Eigen::Matrix4d> poses = depth_parse::GetDepthPoses("data/depth/poses");

    ASSERT_EQ(pcs.size(), 28);    ASSERT_EQ(pcs.size(), poses.size());


    // Poses equal to Python result
    // PCs equal to Python result

    std::vector<cilantro::PointCloud3d> pcs_CLOUDS(pcs.size());
    for (size_t i = 0; i < pcs.size(); ++i)
        pcs_CLOUDS[i] = cilantro::PointCloud3d(pcs[i]);

    auto pc_map = metrics_utils::AggregateMap({pcs_CLOUDS.begin(), pcs_CLOUDS.begin() + 5}, 
                                              {poses.begin(), poses.begin() + 5});

    // std::cout << "MME: " << metrics::GetMME({pcs.begin() + 5, pcs.begin() + 10}, 
    //                                         {poses.begin() + 5, poses.begin() + 10}, 
    //                                         config::CustomConfig(5, 0.2, 30, 5)) 
    //                      << std::endl;


    // std::cout << "MPV: " << metrics::GetMPV(pcs, poses, config::CustomConfig(5, 0.2, 30, 5)) << std::endl;
}
