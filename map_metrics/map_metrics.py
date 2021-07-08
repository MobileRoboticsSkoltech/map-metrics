"""Main module."""


import copy
import math
import networkx as nx
import numpy as np
import open3d as o3d
from sklearn.cluster import AgglomerativeClustering


def aggregate_map(pcs, Ts):
    assert len(pcs) == len(Ts), 'Number of point clouds does not match number of poses'

    Ts = [T @ np.linalg.inv(Ts[0]) for T in Ts]
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


def mpv(pcs, Ts):
    MIN_KNN = 5
    KNN_RAD=1

    pc_map = aggregate_map(pcs, Ts)

    map_tree = o3d.geometry.KDTreeFlann(pc_map)
    points = np.asarray(pc_map.points)

    metric = []
    for i in range(points.shape[0]):
        point = points[i]
        [k, idx, _] = map_tree.search_radius_vector_3d(point, KNN_RAD)
        if len(idx) > MIN_KNN:
            cov = np.cov(points[idx].T)
            eigenvalues = np.linalg.eig(cov)[0]
            metric.append(min(eigenvalues))

    return 0 if len(metric) == 0 else np.mean(metric)


def mom(pcs, Ts, orth_list=None):
    return orth_mme(pcs, Ts, orth_list=None)


def extract_orthogonal_subsets(pc, eps=1e-1, vis=False):
    # Estimate normals if they are not calculated for pc
    if not pc.has_normals():
        pc.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=2, max_nn=30))
    
    pc = pc.uniform_down_sample(4)
    
    normals, lambdas, new_points = build_normals_and_lambdas(pc)

    cut_pcd = o3d.geometry.PointCloud()
    cut_pcd.points = o3d.utility.Vector3dVector(new_points)
    cut_pcd.normals = o3d.utility.Vector3dVector(normals)

    pc = cut_pcd
    # Group normals    
    normals = np.asarray(pc.normals)
    clustering = AgglomerativeClustering(n_clusters=None, linkage="complete", distance_threshold=1e-1
                                     , compute_full_tree=True).fit(normals)
    
    if vis:
        pc_normals = o3d.geometry.PointCloud()
        pc_normals.points = o3d.utility.Vector3dVector(normals)
        o3d.visualization.draw_geometries([pc_normals])

    # Filter out small clusters
    MIN_CLUST_SIZE = 5
    N_clusters = np.unique(clustering.labels_).shape[0]
    labels = clustering.labels_
    huge_clusters = []
    cluster_means = []
    cluster_means_ind = []
    mp = []
    for i in range(N_clusters):
        ind = np.where(labels == i)
        if ind[0].shape[0] > MIN_CLUST_SIZE:
            huge_clusters.append(i)
            cluster_means.append(np.mean(np.vstack(normals)[ind], axis=0))
            cluster_means_ind.append(i)

            if vis:
                pcd = o3d.geometry.PointCloud()
                pcd.points = o3d.utility.Vector3dVector(np.vstack(normals)[ind])
                pcd.paint_uniform_color([0.5 - i / (2 * N_clusters), 0.5 + i / (2 * N_clusters), 
                                         0.5 - i / (2 * N_clusters)])
                mp.append(pcd)
    if vis:
        o3d.visualization.draw_geometries(mp)

    # Normalize means of every cluster
    cluster_means = np.vstack(cluster_means)
    cluster_means = cluster_means / np.linalg.norm(cluster_means, axis=1)[:, None]

    # Generate connectivity graph for normal clusters
    N = cluster_means.shape[0]
    adj_matrix = np.zeros((N, N))
    for i in range(N):
        for j in range(i):
            x = np.abs(np.dot(cluster_means[i], cluster_means[j]))       
            if x < eps:
                adj_matrix[i, j] = 1
                adj_matrix[j, i] = 1
         
    # Find max cliques
    D = nx.Graph(adj_matrix)
    x = nx.algorithms.clique.find_cliques(D)

    # Find cliques with the hugest coverage of points
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
        raise ValueError('AAA')

    max_ind = full_cliques_size.index(max(full_cliques_size))
    max_clique = full_cliques[max_ind]

    # Obtain orth subset and normals for those cliques
    pc_points = np.asarray(pc.points)
    orth_subset = [pc_points[np.where(labels == cluster_means_ind[i])[0]] for i in max_clique]
    pc_normals = np.asarray(pc.normals)
    orth_normals = [pc_normals[np.where(labels == cluster_means_ind[i])[0]] for i in max_clique]
    clique_normals = [cluster_means[i] for i in max_clique]
    return orth_subset, orth_normals, clique_normals


def build_normals_and_lambdas(pc, knn_rad=1):
    pc_tree = o3d.geometry.KDTreeFlann(pc)
    points = np.asarray(pc.points)
    main_normals = np.asarray(pc.normals)
    normals = []
    lambdas = []
    new_points = []
    for i in range(points.shape[0]):
        point = points[i]
        [k, idx, _] = pc_tree.search_radius_vector_3d(point, knn_rad)
        if len(idx) > 3:
            cov = np.cov(points[idx].T)
            eigenvalues, eigenvectors = np.linalg.eig(cov)
            idx = eigenvalues.argsort()
            eigenvalues = eigenvalues[idx]
            eigenvectors = eigenvectors[:, idx]
            if 3 * eigenvalues[0] < eigenvalues[1]:
                normals.append(main_normals[i])
                lambdas.append(eigenvalues[0])
                new_points.append(point)

    return np.vstack(normals), lambdas, np.vstack(new_points)


def orth_mme(pcs, Ts, orth_list=None):
    KNN_RAD = 1
    pc_map = aggregate_map(pcs, Ts)
    map_tree = o3d.geometry.KDTreeFlann(pc_map)
    points = np.asarray(pc_map.points)

    pc = pcs[0]
    pc.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=1.5, max_nn=30))

    if orth_list is None:
        orth_list, _, _ = extract_orthogonal_subsets(pc, eps=1e-1, vis=False)

    pc_map = aggregate_map(pcs, Ts)

    orth_axes_stats = []

    for k, chosen_points in enumerate(orth_list):
        metric = []
        plane_error = []
        for i in range(chosen_points.shape[0]):
            point = chosen_points[i]
            [_, idx, _] = map_tree.search_radius_vector_3d(point, KNN_RAD)
            if len(idx) > 5:
                metric.append(entropy(points[idx]))

        avg_metric = np.median(metric)

        orth_axes_stats.append(avg_metric)

    return np.sum(orth_axes_stats)


def orth_mpv(pcs, Ts, orth_list=None):
    KNN_RAD = 1
    pc_map = aggregate_map(pcs, Ts)
    map_tree = o3d.geometry.KDTreeFlann(pc_map)
    points = np.asarray(pc_map.points)

    pc = pcs[0]
    pc.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=1.5, max_nn=30))

    if orth_list is None:
        orth_list, _, _ = extract_orthogonal_subsets(pc, eps=1e-1, vis=False)

    pc_map = aggregate_map(pcs, Ts)
    orth_axes_stats = []
    
    for k, chosen_points in enumerate(orth_list):
        metric = []
        plane_error = []
        for i in range(chosen_points.shape[0]):
            point = chosen_points[i]
            [_, idx, _] = map_tree.search_radius_vector_3d(point, KNN_RAD)
            if len(idx) > 5:
                metric.append(plane_variance(points[idx]))

        avg_metric = np.median(metric)
    
        orth_axes_stats.append(avg_metric)
    
    return np.sum(orth_axes_stats)


def plane_variance(points):
    cov = np.cov(points.T)
    eigenvalues = np.linalg.eig(cov)[0]
    return min(eigenvalues)


def entropy(points):
    cov = np.cov(points.T)
    det = np.linalg.det(2 * np.pi * np.e * cov)
    if det < 0:
        raise ValueError('Determinant of points is negative!')
    else:
        return 0.5 * np.log(det)
