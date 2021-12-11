# Copyright (c) 2018, Skolkovo Institute of Science and Technology (Skoltech)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#
#  setup.py
#
#  Created on: Jan 22, 2020
#       Author: Lyubov Miloserdova
#               miloslubov@gmail.com
#


import setuptools
import platform, os, ctypes

from setuptools.dist import Distribution
from setuptools import setup, find_packages, Extension


def transform_tag(python, abi, plat):
    if platform.system() == "Darwin":
        python, abi = 'py3', 'none'
        name = plat[:plat.find("_")]
        for i in range(3):
            plat = plat[plat.find("_") + 1:]   # skip name and version of OS
        arch = plat
        version = os.getenv('MACOSX_DEPLOYMENT_TARGET').replace('.', '_')
        plat = name + "_" + version + "_" + arch
    elif platform.system() == "Windows":
        if ctypes.sizeof(ctypes.c_voidp) * 8 > 32:
            plat = "win_" + platform.machine().lower()
        else:
            plat = "win32"
    return python, abi, plat


def wheel_name(**kwargs):
    # create a fake distribution from arguments
    dist = Distribution(attrs=kwargs)
    # finalize bdist_wheel command
    bdist_wheel_cmd = dist.get_command_obj('bdist_wheel')
    bdist_wheel_cmd.ensure_finalized()
    # assemble wheel file name
    distname = bdist_wheel_cmd.wheel_dist_name
    tag = '-'.join(transform_tag(*bdist_wheel_cmd.get_tag()))
    return f'{distname}-{tag}.whl'


setup_kwargs = dict(
    name='map_metrics',
    version='0.0.3',
    packages=find_packages(),
)
file = wheel_name(**setup_kwargs)

setup(**setup_kwargs)
