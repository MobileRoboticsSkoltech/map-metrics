//
// Created by achains on 30.10.2021.
//

#ifndef MAP_METRICS_PC_UTILS_H
#define MAP_METRICS_PC_UTILS_H

#include <vector>

#include "Eigen/Core"


namespace pc_utils{

    Eigen::Matrix4d VectorToEigenMatrix(std::vector<std::vector<double>> vec){
        Eigen::Matrix4d eigen_matrix(4, 4);
        for (Eigen::Index i = 0; i < 4; ++i)
            eigen_matrix.row(i) = Eigen::Map<Eigen::Vector4d>(vec[i].data(), 4);
        return eigen_matrix;
    }

    // TODO: Check return type, should be VectorSet3d
    std::vector<Eigen::Vector3d> VectorToPointCloudPoints(std::vector<double> vec){
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mx(
                vec.data(), static_cast<Eigen::Index>(vec.size() / 4), 4
                );
        Eigen::MatrixX3d mx3d = mx(Eigen::all, Eigen::seq(0, Eigen::last - 1));

        std::vector<Eigen::Vector3d> points(mx3d.rows());
        for (Eigen::Index i = 0; i < mx3d.rows(); ++i)
            points[i] = mx3d.row(i);

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

}

#endif //MAP_METRICS_PC_UTILS_H
