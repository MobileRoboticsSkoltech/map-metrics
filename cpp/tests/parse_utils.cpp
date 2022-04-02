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

    Eigen::Matrix4d VectorToEigenMatrix(std::vector<std::vector<double>> vec){
        Eigen::Matrix4d eigen_matrix(4, 4);
        for (Eigen::Index i = 0; i < 4; ++i)
            eigen_matrix.row(i) = Eigen::Map<Eigen::Vector4d>(vec[i].data(), 4);
        return eigen_matrix;
    }

    cilantro::VectorSet3d VectorToPointCloudPoints(std::vector<double> vec){
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mx(
                vec.data(), static_cast<Eigen::Index>(vec.size() / 4), 4
                );
        Eigen::MatrixX3d mx3d = mx(Eigen::all, Eigen::seq(0, Eigen::last - 1));

        cilantro::VectorSet3d points(3, mx3d.rows());
        for (Eigen::Index i = 0; i < mx3d.rows(); ++i)
            points.col(i) = mx3d.row(i);

        return points;
    }


    std::vector<Eigen::Matrix4d> CalibrateTrajectory(std::vector<std::vector<double>> const & calib,
                                                         std::vector<std::vector<double>> const & trajectory){
        std::vector<Eigen::Matrix4d> calibrated_tj;
        for (auto const & row: trajectory){
            std::vector<std::vector<double>> row_matrix(4, std::vector<double>(4));
            assert (row.size() % 4 == 0);
            for (size_t j = 0; j < row.size(); ++j){
                row_matrix[j / 4][j % 4] = row[j];
            }
            row_matrix[3] = {0, 0, 0, 1};
            Eigen::Matrix4d prod = VectorToEigenMatrix(row_matrix) * VectorToEigenMatrix(calib);
            calibrated_tj.emplace_back(prod);
        }

        return calibrated_tj;
    }

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
}
