===========
Map Metrics for Trajectory Quality
===========

Map metrics toolkit provides a set of metrics **to quantitatively evaluate trajectory quality via estimating 
consistency of the map aggregated from point clouds**.

GPS or Motion Capture systems are not always available in perception systems, or their quality is not enough (GPS on 
small-scale distances) for use as ground truth trajectory. Thus, common full-reference trajectory metrics (APE, 
RPE, and their modifications) could not be applied to evaluate trajectory quality. When 3D sensing technologies (depth 
camera, LiDAR) are available on the perception system, one can alternatively assess trajectory quality --- estimate 
the consistency of the map from registered point clouds via the trajectory.

.. image:: https://img.shields.io/pypi/v/map_metrics.svg
        :target: https://pypi.python.org/pypi/map_metrics

.. image:: https://img.shields.io/travis/anastasiia-kornilova/map_metrics.svg
        :target: https://travis-ci.com/anastasiia-kornilova/map_metrics

.. image:: https://readthedocs.org/projects/map-metrics/badge/?version=latest
        :target: https://map-metrics.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status


.. image:: https://pyup.io/repos/github/anastasiia-kornilova/map_metrics/shield.svg
     :target: https://pyup.io/repos/github/anastasiia-kornilova/map_metrics/
     :alt: Updates

Features
--------
Our toolkit provides implementation of the next metrics:

* Mean Map Entropy (MME), Mean Plane Variance(MPV) [#]_ [#]_
* Mutually Orthogonal Metric (MOM) [#]_ -- has strong correlation with RPE


Documentation: https://map-metrics.readthedocs.io.

Citation
--------

If you use this toolkit or MOM-metric results, please, cite our work:

.. code-block::

    @misc{kornilova2021benchmark,
        title={Be your own Benchmark: No-Reference Trajectory Metric on Registered Point Clouds}, 
        author={Anastasiia Kornilova and Gonzalo Ferrer},
        year={2021},
        eprint={2106.11351},
        archivePrefix={arXiv},
        primaryClass={cs.RO}
    }


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

Links
-----

.. [#] Droeschel, David, Jörg Stückler, and Sven Behnke. "Local multi-resolution representation for 6D motion estimation and mapping with a continuously rotating 3D laser scanner." 2014 IEEE International Conference on Robotics and Automation (ICRA). IEEE, 2014.
.. [#] Razlaw, Jan, et al. "Evaluation of registration methods for sparse 3D laser scans." 2015 European Conference on Mobile Robots (ECMR). IEEE, 2015. 
.. [#] Kornilova, Anastasiia, and Gonzalo Ferrer. "Be your own Benchmark: No-Reference Trajectory Metric on Registered Point Clouds." arXiv preprint arXiv:2106.11351 (2021).