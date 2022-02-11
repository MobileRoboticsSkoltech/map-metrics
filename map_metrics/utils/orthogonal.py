import numpy as np
import networkx as nx
import open3d as o3d

from sklearn.cluster import AgglomerativeClustering


__all__ = ["extract_orthogonal_subsets"]


# TODO: knn_rad <- Config
def _build_normals_and_lambdas(pc, knn_rad=0.2):
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


# TODO: radius, max_nn <- Config
def _estimate_normals(pc):
    if not pc.has_normals():
        pc.estimate_normals(
            search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=0.2, max_nn=40)
        )

    normals, lambdas, new_points = _build_normals_and_lambdas(pc)

    cut_pcd = o3d.geometry.PointCloud()
    cut_pcd.points = o3d.utility.Vector3dVector(new_points)
    cut_pcd.normals = o3d.utility.Vector3dVector(normals)

    return cut_pcd


# TODO: Remove vis
def _filter_clusters(clustering, normals, min_clust_size=5, vis=None):
    if vis is None:
        vis = {}
    n_clusters = np.unique(clustering.labels_).shape[0]
    labels = clustering.labels_
    huge_clusters = []
    cluster_means, cluster_means_ind = [], []
    mp = []

    for i in range(n_clusters):
        ind = np.where(labels == i)
        if ind[0].shape[0] > min_clust_size:
            huge_clusters.append(i)
            cluster_means.append(np.mean(np.vstack(normals)[ind], axis=0))
            cluster_means_ind.append(i)

            if vis:
                pcd = o3d.geometry.PointCloud()
                pcd.points = o3d.utility.Vector3dVector(np.vstack(normals)[ind])
                pcd.paint_uniform_color(
                    [
                        0.5 - i / (2 * n_clusters),
                        0.5 + i / (2 * n_clusters),
                        0.5 - i / (2 * n_clusters),
                    ]
                )
                mp.append(pcd)
    if vis.get("Vis_MP"):
        o3d.visualization.draw_geometries(mp)

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


def extract_orthogonal_subsets(pc, eps=1e-1, vis=None):
    if vis is None:
        vis = {}
    cut_pc = _estimate_normals(pc)
    if vis.get("Vis_Cut_Pc"):
        o3d.visualization.draw_geometries([cut_pc])

    normals = np.asarray(cut_pc.normals)
    # if vis.get("Vis_Normals"):
    #     visualize_normals(normals)

    clustering = AgglomerativeClustering(
        n_clusters=None, distance_threshold=1e-1, compute_full_tree=True
    ).fit(normals)
    labels = clustering.labels_

    cluster_means, cluster_means_ind = _filter_clusters(clustering, normals, vis=vis)

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
