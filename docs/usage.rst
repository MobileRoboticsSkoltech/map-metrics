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

To use it, provide `pcs` --- a list of point clouds pcs in form of `open3d.PointCloud` and `Ts` --- a list of 
corresponding poses in the trajectory in form of `4x4` transformation matrices.

.. code-block:: python
    
    map_metrics.mme(pcs, Ts)


Mean Plane Variance (MPV)
----------------------

Mean Plane Variance calculcates average plane variance of every point vicinity in the aggregated map:

.. math::
    V(P) = \frac{1}{|P|}\sum_{k=1}^{|P|}v(p_k) = \frac{1}{|P|}\sum_{k=1}^{|P|} \lambda_{min}

To use it, provide `pcs` --- a list of point clouds pcs in form of `open3d.PointCloud` and `Ts` --- a list of 
corresponding poses in the trajectory in form of `4x4` transformation matrices.

.. code-block:: python
    
    map_metrics.mpv(pcs, Ts)
