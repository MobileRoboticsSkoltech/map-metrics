from pathlib import Path

from map_metrics.config import LidarConfig
from map_metrics.metrics import mme, mpv, mom

import csv
import numpy as np
import open3d as o3d
import os

import pytest


@pytest.fixture
def calibration_mx():
    return np.array(
        [
            4.276802385584e-04,
            -9.999672484946e-01,
            -8.084491683471e-03,
            -1.198459927713e-02,
            -7.210626507497e-03,
            8.081198471645e-03,
            -9.999413164504e-01,
            -5.403984729748e-02,
            9.999738645903e-01,
            4.859485810390e-04,
            -7.206933692422e-03,
            -2.921968648686e-01,
            0,
            0,
            0,
            1,
        ]
    ).reshape(4, 4)


@pytest.fixture
def poses_path():
    return Path("data/lidar/poses.txt")


@pytest.fixture
def pcs_path():
    return Path("data/lidar/pcs")


@pytest.fixture
def trajectories(poses_path, calibration_mx):
    ts_gt = []
    with open(poses_path, newline="") as csvfile:
        odometry_reader = csv.reader(csvfile, delimiter=" ")
        for row in odometry_reader:
            row = [float(i) for i in row]
            t = np.eye(4)
            t[:3, :4] = np.array(row).reshape(3, 4)
            ts_gt.append(t @ calibration_mx)

    return ts_gt


@pytest.fixture
def pointclouds(pcs_path):
    pc_names = os.listdir(pcs_path)
    pc_names.sort()

    pcs = []
    for _, pc_name in enumerate(pc_names):
        pc = o3d.geometry.PointCloud()
        pc.points = o3d.utility.Vector3dVector(
            np.fromfile(os.path.join(pcs_path, pc_name), dtype=np.float32).reshape(
                -1, 4
            )[:, :3]
        )
        pcs.append(pc.voxel_down_sample(voxel_size=0.05))

    return pcs


# TODO: FINISH TESTS
def test_true():
    assert True


# @pytest.mark.parametrize(
#     "config, metric, expected",
#     [
#         # (LidarConfig, mme, 0.0),
#         # (LidarConfig, mpv, 0.0),
#         # (LidarConfig, mom, 0.0)
#     ],
# )
# def test_metrics(pointclouds, trajectories, config, metric, expected):
#     actual_result = metric(pcs=pointclouds[:3], ts=trajectories[:3], config=config)
#     print(actual_result)
#     assert True
