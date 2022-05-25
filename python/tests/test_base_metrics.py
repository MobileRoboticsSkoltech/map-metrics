import open3d as o3d
import numpy as np
import os

import map_metrics


import pytest


@pytest.fixture
def poses():
    path = 'python/tests/data/depth/poses'
    poses_names = sorted(os.listdir(path))
    return [np.loadtxt(f"{path}/{name}", usecols=range(4)) for name in poses_names]


@pytest.fixture
def points():
    def pcs_to_point_set(pcs):
        return list(map(lambda pc: np.asarray(pc.points).T, pcs))
    path = 'python/tests/data/depth/pcs'
    pcs_names = sorted(os.listdir(path))
    pcs = [o3d.io.read_point_cloud(f"{path}/{name}") for name in pcs_names]
    return pcs_to_point_set(pcs)


@pytest.fixture
def cfg():
    return map_metrics.config.CustomConfig(5, 0.2, 30, 5)


def test_mpv(points, poses, cfg):
    result = map_metrics.mpv(points, poses, cfg)
    assert abs(result - 0.0025390347963358067) < 1e8


def test_mme(points, poses, cfg):
    result = map_metrics.mme(points, poses, cfg)
    assert abs(result - (-3.6144387057257523)) < 1e8
