import numpy as np
import os

# import map_metrics
from map_metrics import map_metrics

import pytest


@pytest.fixture
def poses():
    path = 'python/tests/data/depth/poses'
    poses_names = sorted(os.listdir(path))
    return [np.loadtxt(f"{path}/{name}", usecols=range(4)) for name in poses_names]


@pytest.fixture
def points():
    path = 'python/tests/data/depth/points'
    points_names = sorted(os.listdir(path))
    return [np.load(f"{path}/{name}") for name in points_names]


@pytest.fixture
def cfg():
    return map_metrics.config.CustomConfig(5, 0.2, 30, 5)


def test_mpv(points, poses, cfg):
    result = map_metrics.mpv(points, poses, cfg)
    assert abs(result - 0.0025390347963358067) < 1e8


def test_mme(points, poses, cfg):
    result = map_metrics.mme(points, poses, cfg)
    assert abs(result - (-3.6144387057257523)) < 1e8
