#include "gtest/gtest.h"

#include <iostream>
#include <vector>
#include <Eigen/Core>

#include <cilantro/utilities/point_cloud.hpp>

int main(int argc, char* argv[]){
    testing::InitGoogleTest(&argc, argv);

//    std::vector<Eigen::Vector3d> points_vector(3, Eigen::Vector3d(1, 2, 3));
//
//    std::cout << points_vector[0] << std::endl;
//
//    cilantro::VectorSet3d points_vector_set(3, 4);
//
//    points_vector_set << 1.,  3.,   5.,   2.,
//                         4.,  6.,   1.,   5.,
//                         76., 131., 131., 11;
//
//    std::cout << "\n=======\n" << points_vector_set.row(0) << "\n=======\n";
//    std::cout << "\n=======\n" << points_vector_set.col(0) << "\n=======\n";

    return RUN_ALL_TESTS();
}
