"""Main module."""


import copy
from open3d import o3d


def aggregate_map(pcs, Ts):
    assert len(pcs) == len(Ts), 'Number of point clouds does not match number of poses'

    pc_map = o3d.geometry.PointCloud()
    for i, pc in enumerate(pcs):
        pc_map += copy.deepcopy(pc).transform(Ts[i])

    return pc_map


def mme(pcs, Ts):
    MIN_KNN = 5
    KNN_RAD = 1

    pc_map = aggregate_map(pcs, Ts)

    map_tree = o3d.geometry.KDTreeFlann(pc_map)
    points = np.asarray(pc_map.points)
    metric = []
    for i in range(points.shape[0]):
        point = points[i]
        [k, idx, _] = map_tree.search_radius_vector_3d(point, KNN_RAD)
        if len(idx) > MIN_KNN:
            cov = np.cov(points[idx].T)
            det = np.linalg.det(2 * np.pi * np.e * cov)
            if det > 0:
                metric.append(0.5 * np.log(det))

    return 0 if len(metric) == 0 else np.mean(metric)
