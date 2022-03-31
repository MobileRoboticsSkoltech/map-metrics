#ifndef MAP_METRICS_PARSE_UTILS_H
#define MAP_METRICS_PARSE_UTILS_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <filesystem>

#include <cilantro/core/kd_tree.hpp>
#include <cilantro/utilities/point_cloud.hpp>

namespace parse_utils{
    using PointCloud = cilantro::PointCloud3d;
    using KDTree = cilantro::KDTree3d<>;

    void WritePCToFile(char const * filename, PointCloud const & PC);

    void WriteTjToFile(char const * filename, std::vector<Eigen::Matrix4d> const & Tj);

    std::vector<double> SplitString(std::string const & line, char delim);

    std::vector<std::vector<double>> ParseTrajectory(std::filesystem::path const & path);

    std::vector<double> ParseCSV(std::filesystem::path const & path);

    std::vector<std::vector<double>> GetPointClouds(std::filesystem::path const & path);

    Eigen::Matrix4d VectorToEigenMatrix(std::vector<std::vector<double>> vec);

    cilantro::VectorSet3d VectorToPointCloudPoints(std::vector<double> vec);

    std::vector<Eigen::Matrix4d> CalibrateTrajectory(std::vector<std::vector<double>> const & calib,
                                                     std::vector<std::vector<double>> const & trajectory);
}

#endif //MAP_METRICS_PARSE_UTILS_H