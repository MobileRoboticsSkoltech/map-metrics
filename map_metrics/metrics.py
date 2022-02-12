import copy
import numpy as np
import open3d as o3d

from typing import Optional, Type, Any, List
from nptyping import NDArray

from utils.orthogonal import extract_orthogonal_subsets
from config import BaseConfig, LidarConfig

__all__ = ["aggregate_map", "mme", "mpv", "mom"]


def aggregate_map(pcs: List[o3d.geometry.PointCloud], ts: List[NDArray[(4, 4), np.float64]]) -> o3d.geometry.PointCloud:
    """
    Build a map from point clouds with their poses

    Parameters
    ----------
    pcs: List[o3d.geometry.PointCloud]
        Point Clouds obtained from sensors
    ts: List[NDArray[(4, 4), np.float64]]
        Transformation matrices list (i.e., Point Cloud poses)

    Returns
    -------
    pc_map: o3d.geometry.PointCloud
        Map aggregated from point clouds

    Raises
    ------
    ValueError
        If number of point clouds does not match number of poses
    """
    if len(pcs) != len(ts):
        raise ValueError("Number of point clouds does not match number of poses")

    ts = [np.linalg.inv(ts[0]) @ T for T in ts]
    pc_map = o3d.geometry.PointCloud()
    for i, pc in enumerate(pcs):
        pc_map += copy.deepcopy(pc).transform(ts[i])

    return pc_map


def _plane_variance(points: NDArray[(Any, 3), np.float64]) -> float:
    """
    Compute plane variance of given points

    Parameters
    ----------
    points: NDArray[(Any, 3), np.float64]
        Point Cloud points

    Returns
    -------
    plane_variance: float
        Points plane variance
    """
    cov = np.cov(points.T)
    eigenvalues = np.linalg.eig(cov)[0]
    return min(eigenvalues)


def _entropy(points: NDArray[(Any, 3), np.float64]) -> Optional[float]:
    """
    Compute entropy of given points

    Parameters
    ----------
    points: NDArray[(Any, 3), np.float64]
        Point Cloud points

    Returns
    -------
    entropy: Optional[float]
        Points entropy
    """
    cov = np.cov(points.T)
    det = np.linalg.det(2 * np.pi * np.e * cov)
    if det > 0:
        return 0.5 * np.log(det)

    return None


def _mean_map_metric(pcs, ts, min_knn=5, knn_rad=1, alg=_plane_variance) -> float:
    pc_map = aggregate_map(pcs, ts)

    map_tree = o3d.geometry.KDTreeFlann(pc_map)
    points = np.asarray(pc_map.points)
    metric = []
    for i in range(points.shape[0]):
        point = points[i]
        _, idx, _ = map_tree.search_radius_vector_3d(point, knn_rad)
        if len(idx) > min_knn:
            metric_value = alg(points[idx])
            if metric_value is not None:
                metric.append(metric_value)

    return 0.0 if len(metric) == 0 else np.mean(metric)


def _orth_mpv(pcs, ts, min_knn=5, knn_rad=1, orth_list=None):
    pc_map = aggregate_map(pcs, ts)
    map_tree = o3d.geometry.KDTreeFlann(pc_map)
    points = np.asarray(pc_map.points)

    pc = pcs[0]
    pc.estimate_normals(
        search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=1.5, max_nn=30)
    )

    if orth_list is None:
        orth_list, _, _ = extract_orthogonal_subsets(pc, eps=1e-1, vis=False)

    orth_axes_stats = []

    for k, chosen_points in enumerate(orth_list):
        metric = []
        for i in range(chosen_points.shape[0]):
            point = chosen_points[i]
            [_, idx, _] = map_tree.search_radius_vector_3d(point, knn_rad)
            if len(idx) > min_knn:
                metric.append(_plane_variance(points[idx]))

        avg_metric = np.median(metric)
        orth_axes_stats.append(avg_metric)

    return np.sum(orth_axes_stats)


def mme(pcs, ts, config: Type[BaseConfig] = LidarConfig) -> float:
    return _mean_map_metric(pcs, ts, config.MIN_KNN, config.KNN_RAD, alg=_entropy)


def mpv(pcs, ts, config: Type[BaseConfig] = LidarConfig) -> float:
    return _mean_map_metric(pcs, ts, config.MIN_KNN, config.KNN_RAD, alg=_plane_variance)


def mom(pcs, ts, orth_list=None, config: Type[BaseConfig] = LidarConfig):
    return _orth_mpv(pcs, ts, config.MIN_KNN, config.KNN_RAD, orth_list=orth_list)
