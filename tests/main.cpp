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


//    cilantro::RigidTransform3f tf_ref;
//    Eigen::Matrix3f tmp;
//    tmp = Eigen::AngleAxisf(-0.1f ,Eigen::Vector3f::UnitZ()) *
//          Eigen::AngleAxisf(0.1f, Eigen::Vector3f::UnitY()) *
//          Eigen::AngleAxisf(-0.1f, Eigen::Vector3f::UnitX());
//    tf_ref.linear() = tmp;
//    tf_ref.translation() = Eigen::Vector3f(-0.20f, -0.05f, 0.10f);
//
//    std::cout << std::endl << tf_ref.matrix() << std::endl;


//    cilantro::VectorSet3d points_init(3, 4);
//
//    points_init <<   1.,  3.,   5.,   2.,
//                     4.,  6.,   1.,   5.,
//                     76., 131., 131., 11;
//
//    cilantro::VectorSet3d points_to_add(3, 1);
//
//    points_to_add << 9.,
//                     8.,
//                     7.;
//    Eigen::Index old_size = points_init.cols();
//    points_init.conservativeResize(3, 5);
//    for (Eigen::Index i = 0; i < points_to_add.cols(); ++i) points_init.col(old_size + i) = points_to_add.col(i);
//
////    cilantro::VectorSet3d new_points(3, points_init.cols() + points_to_add.cols());
////    for (Eigen::Index i = 0; i < points_init.cols(); ++i) new_points.col(i) = points_init.col(i);
////    for (Eigen::Index i = 0; i < points_to_add.cols(); ++i) new_points.col(points_init.cols() + i) = points_to_add.col(i);
//
//    std::cout << points_init << std::endl;

    return RUN_ALL_TESTS();
}
