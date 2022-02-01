=====
Usage
=====

To use Map Metrics in a project::

    import map_metrics


Mean Map Entropy (MME)
----------------------

Mean Map Entropy calculcates average entropy of every point vicinity in the aggregated map:

.. math::

    h(p_k) = \frac{1}{2}\det \big(2\pi e \Sigma(W(p_k)) \big)

    H(P) = \frac{1}{|P|}\sum_{k=1}^{|P|} h(q_k)

To use it, provide `pcs` --- a list of point clouds pcs in the form of `open3d.PointCloud` and `Ts` --- a list of
corresponding poses in the trajectory in the form of `4x4` transformation matrices.

.. code-block:: python

    map_metrics.mme(pcs, Ts)


Mean Plane Variance (MPV)
-------------------------

Mean Plane Variance calculcates average plane variance of every point vicinity in the aggregated map:

.. math::
    V(P) = \frac{1}{|P|}\sum_{k=1}^{|P|}v(p_k) = \frac{1}{|P|}\sum_{k=1}^{|P|} \lambda_{min}

To use it, provide `pcs` --- a list of point clouds pcs in the form of `open3d.PointCloud` and `Ts` --- a list of
corresponding poses in the trajectory in the form of `4x4` transformation matrices.

.. code-block:: python

    map_metrics.mpv(pcs, Ts)


Mutually Orthogonal Metric (MOM)
--------------------------------

Mutual Orthogonality is a concept of considering not all points in the map but only ones from mutually orthogonal
surfaces. Mean Plane Variance over those points provides (as described in our paper) correlation with Relative Pose
Error (RPE) --- one of the popular full-reference metrics for trajectories.

To use it, provide `pcs` --- a list of point clouds pcs in the form of `open3d.PointCloud` and `Ts` --- a list of
corresponding poses in the trajectory in the form of `4x4` transformation matrices.

.. code-block:: python

    map_metrics.mom(pcs, Ts)

The default usage of the method assumes extraction of points from mutually orthogonal surfaces during the method
execution and therefore increase calculation time. One can extract those points manually one time for specific set
of point cloud and use it as parameter to calculate MOM faster.

.. code-block:: python

    orth_list, _, _ = map_metrics.extract_orthogonal_subsets(pcs[0])
    print(map_metrics.orth_mme(pcs, Ts_gt, orth_list=orth_list))
