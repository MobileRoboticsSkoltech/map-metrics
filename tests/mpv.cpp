#include <gtest/gtest.h>

#include "test_utils.h"

TEST(Parser, ParsePointClouds){
    ASSERT_EQ(test_utils::GetPointClouds("data/kitti_00").size(), 20);
}