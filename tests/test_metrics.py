from pathlib import Path

from map_metrics.config import DepthConfig
from map_metrics.metrics import mme, mpv, mom
from map_metrics.utils.orthogonal import read_orthogonal_subset

import numpy as np
import open3d as o3d
import os

import pytest


@pytest.fixture
def depth_trajectories():
    ts_folder = Path("tests/data/depth/poses")
    ts_names = sorted(os.listdir(ts_folder))
    trajectories = []
    for name in ts_names:
        trajectories.append(np.loadtxt(f"{ts_folder}/{name}", usecols=range(4)))

    return trajectories


@pytest.fixture
def depth_pointclouds():
    pcs_folder = Path("tests/data/depth/pcs")
    pc_names = sorted(os.listdir(pcs_folder))
    pcs = []
    for name in pc_names:
        pcs.append(o3d.io.read_point_cloud(f"{pcs_folder}/{name}"))

    return pcs


@pytest.fixture
def depth_orthsubset(depth_trajectories):
    orth_list_name = "tests/data/depth/orth_subset/orth-0091.npy"
    orth_tj_name = "tests/data/depth/poses/pose-0091.txt"

    orth_list = read_orthogonal_subset(
        orth_subset_name=Path(orth_list_name),
        orth_pose_name=Path(orth_tj_name),
        ts=depth_trajectories,
    )

    return orth_list


@pytest.mark.parametrize(
    "config, metric, expected",
    [
        (DepthConfig, mme, -3.478778715),
        (DepthConfig, mpv, 0.003242216),
    ],
)
def test_basic_metrics(depth_pointclouds, depth_trajectories, config, metric, expected):
    actual_result = metric(pcs=depth_pointclouds, ts=depth_trajectories, config=config)
    assert abs(actual_result - expected) < 1e-3


def test_mom(depth_pointclouds, depth_trajectories, depth_orthsubset):
    actual_result = mom(
        depth_pointclouds, depth_trajectories, depth_orthsubset, config=DepthConfig
    )
    assert abs(actual_result - 0.007381899) < 1e-3
