import copy
import numpy as np
import open3d as o3d

from typing import Optional
from utils.orthogonal import extract_orthogonal_subsets


__all__ = ["aggregate_map", "mme", "mpv", "mom"]


def aggregate_map(pcs, ts):
    assert len(pcs) == len(ts), "Number of point clouds does not match number of poses"

    ts = [np.linalg.inv(ts[0]) @ T for T in ts]
    pc_map = o3d.geometry.PointCloud()
    for i, pc in enumerate(pcs):
        pc_map += copy.deepcopy(pc).transform(ts[i])

    return pc_map


def _plane_variance(points) -> float:
    cov = np.cov(points.T)
    eigenvalues = np.linalg.eig(cov)[0]
    return min(eigenvalues)


def _entropy(points) -> Optional[float]:
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
    # TODO: radius, max_nn -- from config
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


# TODO: min_knn, knn_rad <- config
def mme(pcs, ts, min_knn=5, knn_rad=1) -> float:
    return _mean_map_metric(pcs, ts, min_knn, knn_rad, alg=_entropy)


# TODO: min_knn, knn_rad <- config
def mpv(pcs, ts, min_knn=5, knn_rad=1) -> float:
    return _mean_map_metric(pcs, ts, min_knn, knn_rad, alg=_plane_variance)


def mom(pcs, ts, orth_list=None):
    return _orth_mpv(pcs, ts, orth_list=orth_list)
