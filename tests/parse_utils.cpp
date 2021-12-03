//
// Created by achains on 30.10.2021.
//

#include <algorithm>
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

//    void WritePCToFile(char const * filename, PointCloud const & PC){
//        std::ofstream output(filename);
//        for (Eigen::Vector3d const & row : PC.points_){
//            for (int i = 0; i < 3; ++i) output << std::setprecision(17) << row[i] << ' ';
//            output << '\n';
//        }
//    }
//
//    void WriteTjToFile(char const * filename, std::vector<Eigen::Matrix4d> const & Tj){
//        std::ofstream output(filename);
//        for (Eigen::Matrix4d const & mx : Tj){
//            for (int i = 0; i < 4; ++i){
//                for (int j = 0; j < 4; ++j) output << std::setprecision(17) << mx(i, j) << ' ';
//                output << '\n';
//            }
//        }
//    }

    std::vector<double> SplitString(std::string const & line, char delim){
        std::istringstream ss(line);
        std::string token;
        std::vector<double> numbers;
        while(std::getline(ss, token, delim)){
            numbers.push_back(std::stod(token));
        }

        return numbers;
    }

    std::vector<std::vector<double>> ParseTrajectory(std::filesystem::path const & path){
        std::ifstream data(path);
        std::vector<std::vector<double>> trajectory;
        std::string line;

        while(std::getline(data, line)){
            trajectory.push_back(SplitString(line, ' '));
        }

        return trajectory;
    }

    std::vector<double> ParseCSV(std::filesystem::path const & path){
        std::ifstream data(path, std::ios::binary);
        std::vector<double> PC;
        char* line = new char[sizeof(float)];
        while(data.read(line, sizeof(float))){
            PC.push_back(*reinterpret_cast<float*>(line));
        }
        delete[] line;

        return PC;
    }

    std::vector<std::vector<double>> GetPointClouds(std::filesystem::path const & path){
        std::vector<std::vector<double>> point_clouds;
        std::vector<std::filesystem::path> PC_file_names;
        for (const auto & entry : std::filesystem::directory_iterator(path)){
            PC_file_names.push_back(entry.path());
        }
        std::sort(PC_file_names.begin(), PC_file_names.end());
        for (const auto & filepath : PC_file_names){
            point_clouds.emplace_back(ParseCSV(filepath));
        }

        return point_clouds;
    }
}
