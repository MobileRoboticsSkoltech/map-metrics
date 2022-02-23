import copy
import numpy as np
import open3d as o3d

from typing import Optional, Type, Any, List, Callable
from nptyping import NDArray

from map_metrics.utils.orthogonal import extract_orthogonal_subsets
from map_metrics.config import BaseConfig, LidarConfig

__all__ = ["aggregate_map", "mme", "mpv", "mom"]


def aggregate_map(
    pcs: List[o3d.geometry.PointCloud], ts: List[NDArray[(4, 4), np.float64]]
) -> o3d.geometry.PointCloud:
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


def _mean_map_metric(
    pcs: List[o3d.geometry.PointCloud],
    ts: List[NDArray[(4, 4), np.float64]],
    config: Type[BaseConfig] = LidarConfig,
    alg: Callable = _plane_variance,
) -> float:
    """
    No-reference metric algorithms helper

    Parameters
    ----------
    pcs: List[o3d.geometry.PointCloud]
        Point Clouds obtained from sensors
    ts: List[NDArray[(4, 4), np.float64]]
        Transformation matrices list (i.e., Point Cloud poses)
    config: BaseConfig
        Scene hyperparameters
    alg: Callable
        Metric algorithm basis (e.g., plane variance, entropy)
    Returns
    -------
    mean: float
        Mean of given metric algorithm values
    """
    pc_map = aggregate_map(pcs, ts)

    map_tree = o3d.geometry.KDTreeFlann(pc_map)
    points = np.asarray(pc_map.points)
    metric = []
    for i in range(points.shape[0]):
        point = points[i]
        _, idx, _ = map_tree.search_radius_vector_3d(point, config.KNN_RAD)
        if len(idx) > config.MIN_KNN:
            metric_value = alg(points[idx])
            if metric_value is not None:
                metric.append(metric_value)

    return 0.0 if len(metric) == 0 else np.mean(metric)


def _orth_mpv(
    pcs: List[o3d.geometry.PointCloud],
    ts: List[NDArray[(4, 4), np.float64]],
    config: Type[BaseConfig] = LidarConfig,
    orth_list: List[o3d.geometry.PointCloud] = None,
):
    """

    Parameters
    ----------
    pcs: List[o3d.geometry.PointCloud]
        Point Clouds obtained from sensors
    ts: List[NDArray[(4, 4), np.float64]]
        Transformation matrices list (i.e., Point Cloud poses)
    config: BaseConfig
        Scene hyperparameters
    orth_list: List[o3d.geometry.PointCloud], default=None
        List of orthogonal planes of the map

    Returns
    -------
    val: float
        The value of MPV computed on orthogonal planes of the map
    """
    pc_map = aggregate_map(pcs, ts)
    map_tree = o3d.geometry.KDTreeFlann(pc_map)
    points = np.asarray(pc_map.points)

    if orth_list is None:
        pc = pcs[0]
        orth_list, _, _ = extract_orthogonal_subsets(pc, config=config, eps=1e-1)

    orth_axes_stats = []

    for chosen_points in orth_list:
        metric = []
        for i in range(np.asarray(chosen_points).shape[0]):
            point = chosen_points[i]
            _, idx, _ = map_tree.search_radius_vector_3d(point, config.KNN_RAD)
            # TODO: add 3 to config
            if len(idx) > 3:
                metric.append(_plane_variance(points[idx]))

        avg_metric = np.median(metric)
        orth_axes_stats.append(avg_metric)

    return np.sum(orth_axes_stats)


def mme(
    pcs: List[o3d.geometry.PointCloud],
    ts: List[NDArray[(4, 4), np.float64]],
    config: Type[BaseConfig] = LidarConfig,
) -> float:
    """
    Mean Map Entropy
    A no-reference metric algorithm based on entropy

    Parameters
    ----------
    pcs: List[o3d.geometry.PointCloud]
        Point Clouds obtained from sensors
    ts: List[NDArray[(4, 4), np.float64]]
        Transformation matrices list (i.e., Point Cloud poses)
    config: BaseConfig
        Scene hyperparameters

    Returns
    -------
    mean: float
        Mean of given metric algorithm values
    """
    return _mean_map_metric(pcs, ts, config, alg=_entropy)


def mpv(
    pcs: List[o3d.geometry.PointCloud],
    ts: List[NDArray[(4, 4), np.float64]],
    config: Type[BaseConfig] = LidarConfig,
) -> float:
    """
    Mean Plane Variance
    A no-reference metric algorithm based on plane variance

    Parameters
    ----------
    pcs: List[o3d.geometry.PointCloud]
        Point Clouds obtained from sensors
    ts: List[NDArray[(4, 4), np.float64]]
        Transformation matrices list (i.e., Point Cloud poses)
    config: BaseConfig
        Scene hyperparameters

    Returns
    -------
    mean: float
        Mean of given metric algorithm values
    """
    return _mean_map_metric(pcs, ts, config, alg=_plane_variance)


def mom(
    pcs: List[o3d.geometry.PointCloud],
    ts: List[NDArray[(4, 4), np.float64]],
    orth_list: List[o3d.geometry.PointCloud] = None,
    config: Type[BaseConfig] = LidarConfig,
):
    """
    Mutually Orthogonal Metric
    A no-reference metric algorithm based on MPV on orthogonal planes subset

    Parameters
    ----------
    pcs: List[o3d.geometry.PointCloud]
        Point Clouds obtained from sensors
    ts: List[NDArray[(4, 4), np.float64]]
        Transformation matrices list (i.e., Point Cloud poses)
    orth_list: List[o3d.geometry.PointCloud], default=None
        List of orthogonal planes of the map
    config: BaseConfig
        Scene hyperparameters

    Returns
    -------
    mean: float
        Mean of given metric algorithm values
    """
    return _orth_mpv(pcs, ts, config, orth_list=orth_list)
