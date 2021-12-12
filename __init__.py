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
#  __init__.py
#
#  Created on: Nov 11, 2020
#       Author: Lyubov Miloserdova
#               miloslubov@gmail.com
#

try:
    from map_metrics import map_metrics

    for module in dir(map_metrics):
        n = len(module) - 1
        if not (module[:2] == '__' and module[n:n-2:-1] == '__') and module.count('.') == 0:
            globals()[module] = getattr(map_metrics, module)

    del map_metrics
except ImportError:
    import platform

    if platform.system() == "Windows":
        import subprocess
        import ctypes

        if ctypes.sizeof(ctypes.c_voidp) * 8 > 32:
            msvc_not_instaled = min(
                subprocess.call(
                    'REG QUERY "HKEY_LOCAL_MACHINE\\SOFTWARE\\Wow6432Node\\Microsoft\\VisualStudio\\14.0\\VC\\Runtimes\\x64"',
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                ),
                subprocess.call(
                    'REG QUERY "HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\14.0\\VC\\Runtimes\\x64"',
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                )
            )
        else:
            msvc_not_instaled = min(
                subprocess.call(
                    'REG QUERY "HKEY_LOCAL_MACHINE\\SOFTWARE\\Wow6432Node\\Microsoft\\VisualStudio\\14.0\\VC\\Runtimes\\x86"',
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                ),
                subprocess.call(
                    'REG QUERY "HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\14.0\\VC\\Runtimes\\x86"',
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                )
            )

        if msvc_not_instaled:
            import sys

            sys.tracebacklimit = 0
            raise ImportError(
                "You don't have Microsoft Visual C++ installed. Please follow the link and install redistributable " +
                "package: https://docs.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-160" +
                "#visual-studio-2015-2017-2019-and-2022") from None

    del platform
