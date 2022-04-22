//
// Created by achains on 30.10.2021.
//

#ifndef MAP_METRICS_PC_UTILS_H
#define MAP_METRICS_PC_UTILS_H

#include <vector>
#include "cilantro/core/kd_tree.hpp"

#include "Eigen/Core"


namespace pc_utils{

    Eigen::Matrix4d VectorToEigenMatrix(std::vector<std::vector<double>> vec);

    cilantro::VectorSet3d VectorToPointCloudPoints(std::vector<double> vec);

    std::vector<Eigen::Matrix4d> CalibrateTrajectory(std::vector<std::vector<double>> const & calib,
                                                         std::vector<std::vector<double>> const & trajectory);

}

#endif //MAP_METRICS_PC_UTILS_H
