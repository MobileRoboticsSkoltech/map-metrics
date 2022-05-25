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


try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
    import platform, os, ctypes

    class bdist_wheel(_bdist_wheel):

        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            if platform.system() == "Darwin":
                self.root_is_pure = False

        def get_tag(self):
            python, abi, plat = _bdist_wheel.get_tag(self)
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

except ImportError:
    bdist_wheel = None


setup_kwargs = dict(
    name='map_metrics',
    version='0.0.6',
    packages=find_packages(),
    cmdclass={'bdist_wheel': bdist_wheel}
)

setup(**setup_kwargs)
