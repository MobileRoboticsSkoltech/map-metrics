from pathlib import Path

import csv
import numpy as np
import os

import map_metrics

import pytest


@pytest.fixture
def t_calib():
    return np.array([4.276802385584e-04, -9.999672484946e-01, -8.084491683471e-03, -1.198459927713e-02,
                     -7.210626507497e-03, 8.081198471645e-03, -9.999413164504e-01, -5.403984729748e-02,
                     9.999738645903e-01, 4.859485810390e-04, -7.206933692422e-03, -2.921968648686e-01,
                     0, 0, 0, 1]).reshape(4, 4)


@pytest.fixture
def gt_traj_file():
    return Path('python/tests/data/00.txt')


@pytest.fixture
def pcs_folder():
    return Path('python/tests/data/kitti_00/')


@pytest.fixture
def trajectories(gt_traj_file, t_calib):
    ts_gt = []
    with open(gt_traj_file, newline='') as csvfile:
        odometry_reader = csv.reader(csvfile, delimiter=' ')
        for row in odometry_reader:
            row = [float(i) for i in row]
            t = np.eye(4)
            t[:3, :4] = np.array(row).reshape(3, 4)
            ts_gt.append(t @ t_calib)

    return ts_gt


@pytest.fixture
def pointclouds(pcs_folder):
    pc_names = os.listdir(pcs_folder)
    pc_names.sort()

    pcs = []
    for _, pc_name in enumerate(pc_names):
        points = np.fromfile(os.path.join(pcs_folder, pc_name), dtype=np.float32).reshape(-1, 4)[:, :3].T
        pcs.append(points)


def test_mpv(trajectories, pointclouds):
    result = map_metrics.mpv(trajectories, pointclouds, 5, 1.0)
    assert abs(result - 0.08949105590777887) < 1e8


def test_mme(trajectories, pointclouds):
    result = map_metrics.mme(trajectories, pointclouds, 5, 1.0)
    assert abs(result - 1.2553150544133582) < 1e8
