#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <filesystem>

#include "open3d/Open3D.h"

namespace test_utils{
    using open3d::geometry::PointCloud;
    using open3d::geometry::KDTreeFlann;

    void WritePCToFile(char const * filename, std::vector<double> const & PC){
        std::ofstream output(filename);
        for (double const & number: PC)
            output << std::setprecision(17) << number << '\n';
    }

    std::vector<double> ParseCSV(std::filesystem::path const & path){
        std::ifstream data(path, std::ios::binary);
        std::vector<double> PC;
        char* line = new char[sizeof(double)];
        while(data.read(line, sizeof(double))){
            PC.push_back(*reinterpret_cast<double*>(line));
        }
        delete[] line;

        return PC;
    }

    std::vector<std::vector<double>> GetPointClouds(char const * path){
        std::vector<std::vector<double>> point_clouds;
        for (const auto & entry : std::filesystem::directory_iterator(path)){
            point_clouds.emplace_back(ParseCSV(entry.path()));
        }

        return point_clouds;
    }
}