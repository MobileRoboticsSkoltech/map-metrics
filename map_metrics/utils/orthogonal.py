import numpy as np
import networkx as nx
import open3d as o3d

from typing import Type, List
from nptyping import NDArray

from map_metrics.config import BaseConfig, LidarConfig
from sklearn.cluster import AgglomerativeClustering
from pathlib import Path


__all__ = ["extract_orthogonal_subsets", "read_orthogonal_subset"]


def _build_normals_and_lambdas(pc, knn_rad):
    pc_tree = o3d.geometry.KDTreeFlann(pc)
    points = np.asarray(pc.points)
    main_normals = np.asarray(pc.normals)
    normals = []
    lambdas = []
    new_points = []
    # TODO: Add mutable tqdm bar
    for i in range(points.shape[0]):
        point = points[i]
        _, idx, _ = pc_tree.search_radius_vector_3d(point, knn_rad)
        if len(idx) > 3:
            cov = np.cov(points[idx].T)
            eigenvalues, eigenvectors = np.linalg.eig(cov)
            idx = eigenvalues.argsort()
            eigenvalues = eigenvalues[idx]
            if 100 * eigenvalues[0] < eigenvalues[1]:
                normals.append(main_normals[i])
                lambdas.append(eigenvalues[0])
                new_points.append(point)

    return np.vstack(normals), lambdas, np.vstack(new_points)


def _estimate_normals(pc, knn_rad, max_nn):
    if not pc.has_normals():
        pc.estimate_normals(
            search_param=o3d.geometry.KDTreeSearchParamHybrid(
                radius=knn_rad, max_nn=max_nn
            )
        )

    normals, lambdas, new_points = _build_normals_and_lambdas(pc, knn_rad)

    cut_pcd = o3d.geometry.PointCloud()
    cut_pcd.points = o3d.utility.Vector3dVector(new_points)
    cut_pcd.normals = o3d.utility.Vector3dVector(normals)

    return cut_pcd


def _filter_clusters(clustering, normals, min_clust_size):
    n_clusters = np.unique(clustering.labels_).shape[0]
    labels = clustering.labels_
    huge_clusters = []
    cluster_means, cluster_means_ind = [], []

    for i in range(n_clusters):
        ind = np.where(labels == i)
        if ind[0].shape[0] > min_clust_size:
            huge_clusters.append(i)
            cluster_means.append(np.mean(np.vstack(normals)[ind], axis=0))
            cluster_means_ind.append(i)

    # Normalize means of every cluster
    cluster_means = np.vstack(cluster_means)
    cluster_means = cluster_means / np.linalg.norm(cluster_means, axis=1)[:, None]

    return cluster_means, cluster_means_ind


def _find_max_clique(labels, cluster_means, cluster_means_ind, eps=1e-1):
    N = cluster_means.shape[0]
    adj_matrix = np.zeros((N, N))
    for i in range(N):
        for j in range(i):
            x = np.abs(np.dot(cluster_means[i], cluster_means[j]))
            if x < eps:
                adj_matrix[i, j] = 1
                adj_matrix[j, i] = 1

    D = nx.Graph(adj_matrix)
    x = nx.algorithms.clique.find_cliques(D)

    full_cliques_size = []
    full_cliques = []
    for clique in x:
        if len(clique) > 2:
            amount = 0
            for j in clique:
                amount += np.sum(labels == cluster_means_ind[j])
            full_cliques_size.append(amount)
            full_cliques.append(clique)

    if len(full_cliques) == 0:
        raise ValueError("Length of full_cliques == 0")

    max_ind = full_cliques_size.index(max(full_cliques_size))
    return full_cliques[max_ind]


def extract_orthogonal_subsets(pc, config: Type[BaseConfig] = LidarConfig, eps=1e-1):
    cut_pc = _estimate_normals(pc, knn_rad=config.KNN_RAD, max_nn=config.MAX_NN)

    normals = np.asarray(cut_pc.normals)

    clustering = AgglomerativeClustering(
        n_clusters=None, distance_threshold=1e-1, compute_full_tree=True
    ).fit(normals)
    labels = clustering.labels_

    cluster_means, cluster_means_ind = _filter_clusters(
        clustering, normals, min_clust_size=config.MIN_CLUST_SIZE
    )

    max_clique = _find_max_clique(labels, cluster_means, cluster_means_ind, eps=eps)

    # Obtain orth subset and normals for those cliques
    pc_points = np.asarray(cut_pc.points)
    orth_subset = [
        pc_points[np.where(labels == cluster_means_ind[i])[0]] for i in max_clique
    ]
    pc_normals = np.asarray(cut_pc.normals)
    orth_normals = [
        pc_normals[np.where(labels == cluster_means_ind[i])[0]] for i in max_clique
    ]
    clique_normals = [cluster_means[i] for i in max_clique]

    return orth_subset, orth_normals, clique_normals


def read_orthogonal_subset(
    orth_subset_name: Path, orth_pose_name: Path, ts: List[NDArray[(4, 4), np.float64]]
):
    """Read and aggregate an orthogonal subset

    Parameters
    ----------
    orth_subset_name: Path
        Orthogonal subset data
    orth_pose_name: Path
        Pose of orthogonal subset in the map
    ts: List[NDArray[(4, 4), np.float64]]
        Transformation matrices list (i.e., Point Cloud poses)

    Returns
    -------
    orth_subset
        Aggregated orthogonal subset
    """
    orth_list = np.load(str(orth_subset_name), allow_pickle=True)
    orth_pose = np.loadtxt(str(orth_pose_name), usecols=range(4))
    orth_pose = np.linalg.inv(ts[0]) @ orth_pose

    return [
        o3d.geometry.PointCloud(o3d.utility.Vector3dVector(surface))
        .transform(orth_pose)
        .points
        for surface in orth_list
    ]
