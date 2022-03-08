# Copyright (c) 2022, Skolkovo Institute of Science and Technology (Skoltech)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from abc import ABC


__all__ = ["BaseConfig", "DepthConfig", "LidarConfig", "CustomConfig"]


class BaseConfig(ABC):
    """
    Base Config Class
    """

    KNN_RAD = None
    MIN_KNN = None
    MAX_NN = None
    MIN_CLUST_SIZE = None


class DepthConfig(BaseConfig):
    """
    Config recommended for data obtained from Depth Camera
    """

    KNN_RAD = 0.2
    MIN_KNN = 5
    MAX_NN = 30
    MIN_CLUST_SIZE = 5


class LidarConfig(BaseConfig):
    """
    Config recommended for data obtained from LiDAR
    """

    KNN_RAD = 1
    MIN_KNN = 5
    MAX_NN = 30
    MIN_CLUST_SIZE = 5


class CustomConfig(BaseConfig):
    """
    Configurable Config

    Attributes
    ----------
    KNN_RAD: float, default=1.0
        Estimated radius between a pair of points on the map
        The value is given in meters
    MIN_KNN: int, default=5
        Minimum number of a point neighbors
    MAX_NN: int, default=30
        At most MAX_NN nearest neighbors that have distances to the anchor point less than a given radius
    MIN_CLUST_SIZE: int, default=5
        Minimal acceptable cluster size in orthogonal extraction
    """

    def __init__(self, knn_rad=1.0, min_knn=5, max_nn=30, min_clust_size=5):
        self.KNN_RAD = knn_rad
        self.MIN_KNN = min_knn
        self.MAX_NN = max_nn
        self.MIN_CLUST_SIZE = min_clust_size
