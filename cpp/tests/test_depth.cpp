#include <gtest/gtest.h>
#include <chrono>
#include <iostream>

#include "cilantro/3rd_party/tinyply/tinyply.h"

#include "parse_utils.h"
#include "pc_utils.h"

#include "metrics.h"

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

namespace kdsearch{
    double GetDistance(Eigen::Vector3d const & x, Eigen::Vector3d const & y){
        double distance = sqrt(pow(x[0] - y[0], 2.0) + pow(x[1] - y[1], 2.0) + pow(x[2] - y[2], 2.0));
        return distance;
    }

    std::vector<Eigen::Index> radius_search(cilantro::VectorSet3d const & points, 
                                            Eigen::Vector3d const & query, 
                                            double radius){
        std::vector<Eigen::Index> indices;
        for (Eigen::Index i = 0; i < points.cols(); ++i){
            if (GetDistance(query, points.col(i)) < radius)
                indices.push_back(i);
        }

        return indices;
    }
}

TEST(DepthTest, MME){
    std::vector<cilantro::VectorSet3d> pcs = depth_parse::GetDepthPCs("data/depth/pcs");
    std::vector<Eigen::Matrix4d> poses = depth_parse::GetDepthPoses("data/depth/poses");
    ASSERT_EQ(pcs.size(), 28);    ASSERT_EQ(pcs.size(), poses.size());

    auto start_time = std::chrono::system_clock::now();
    double actual_mme = metrics::GetMME(pcs, poses, config::CustomConfig(5, 0.2, 30, 5));
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);

    std::cout << "MME: " << actual_mme << std::endl;
    std::cout << "Elapsed time (ms): " << elapsed_milliseconds.count() << std::endl;

    ASSERT_LE(std::fabs(actual_mme - -3.6144387057257523), 1e-8);
}

TEST(DepthTest, MPV){
    std::vector<cilantro::VectorSet3d> pcs = depth_parse::GetDepthPCs("data/depth/pcs");
    std::vector<Eigen::Matrix4d> poses = depth_parse::GetDepthPoses("data/depth/poses");
    ASSERT_EQ(pcs.size(), 28);    ASSERT_EQ(pcs.size(), poses.size());

    auto start_time = std::chrono::system_clock::now();
    double actual_mpv = metrics::GetMPV(pcs, poses, config::CustomConfig(5, 0.2, 30, 5));
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);

    std::cout << "MPV: " << actual_mpv << std::endl;
    std::cout << "Elapsed time (ms): " << elapsed_milliseconds.count() << std::endl;

    ASSERT_LE(std::fabs(actual_mpv - 0.0025390347963358565), 1e-8);
}