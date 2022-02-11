from abc import ABC


class Config(ABC):
    pass


class DepthConfig(Config):
    def __init__(self):
        pass


class LidarConfig(Config):
    def __init__(self):
        pass


class CustomConfig(Config):
    def __init__(self):
        pass
