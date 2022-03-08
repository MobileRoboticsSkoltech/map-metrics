# map-metrics
Map metrics toolkit provides a set of metrics **to quantitatively evaluate trajectory quality via estimating
consistency of the map aggregated from point clouds**.

GPS or Motion Capture systems are not always available in perception systems, or their quality is not enough (GPS on
small-scale distances) for use as ground truth trajectory. Thus, common full-reference trajectory metrics (APE,
RPE, and their modifications) could not be applied to evaluate trajectory quality. When 3D sensing technologies (depth
camera, LiDAR) are available on the perception system, one can alternatively assess trajectory quality -- estimate
the consistency of the map from registered point clouds via the trajectory.

## Installation
```bash
$ pip install map-metrics
```

## Usage
Run metric algorithms from `map_metrics` on your point cloud data.

```python
from map_metrics.metrics import mme
from map_metrics.config import LidarConfig

result = mme(pointclouds, poses, config=LidarConfig)
```

## License
License...

## Credits
Credits...

## Citation
If you use this toolkit or MOM-metric results, please, cite our work:

    @misc{kornilova2021benchmark,
        title={Be your own Benchmark: No-Reference Trajectory Metric on Registered Point Clouds},
        author={Anastasiia Kornilova and Gonzalo Ferrer},
        year={2021},
        eprint={2106.11351},
        archivePrefix={arXiv},
        primaryClass={cs.RO}
    }
