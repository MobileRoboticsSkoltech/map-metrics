

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

from map_metrics import map_metrics

for module in dir(map_metrics):
    n = len(module) - 1
    if not (module[:2] == '__' and module[n:n-2:-1] == '__') and module.count('.') == 0:
        globals()[module] = getattr(map_metrics, module)

del map_metrics
