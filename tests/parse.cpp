#include <gtest/gtest.h>

#include "parse_utils.h"

TEST(Parser, ParsePointClouds){
    ASSERT_EQ(parse_utils::GetPointClouds("data/kitti_00").size(), 20);
}

TEST(Parser, ParseTrajectory){
    auto trajectory = parse_utils::ParseTrajectory("data/00.txt");
    ASSERT_EQ(trajectory.size(), 20);
    ASSERT_EQ(trajectory[0].size(), 12);
}