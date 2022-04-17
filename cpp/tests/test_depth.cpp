#include <gtest/gtest.h>
#include <chrono>
#include <iostream>

#include "cilantro/3rd_party/tinyply/tinyply.h"

#include "parse_utils.h"
#include "pc_utils.h"

#include "metrics.h"

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
    std::vector<cilantro::VectorSet3d> pcs = parse_utils::depth_parse::GetDepthPCs("data/depth/pcs");
    std::vector<Eigen::Matrix4d> poses = parse_utils::depth_parse::GetDepthPoses("data/depth/poses");
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
    std::vector<cilantro::VectorSet3d> pcs = parse_utils::depth_parse::GetDepthPCs("data/depth/pcs");
    std::vector<Eigen::Matrix4d> poses = parse_utils::depth_parse::GetDepthPoses("data/depth/poses");
    ASSERT_EQ(pcs.size(), 28);    ASSERT_EQ(pcs.size(), poses.size());

    auto start_time = std::chrono::system_clock::now();
    double actual_mpv = metrics::GetMPV(pcs, poses, config::CustomConfig(5, 0.2, 30, 5));
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);

    std::cout << "MPV: " << actual_mpv << std::endl;
    std::cout << "Elapsed time (ms): " << elapsed_milliseconds.count() << std::endl;

    ASSERT_LE(std::fabs(actual_mpv - 0.0025390347963358565), 1e-8);
}