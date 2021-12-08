#!/bin/bash
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
#  build-wheels-linux.sh
#
#  Created on: Nov 28, 2020
#       Author: Lyubov Miloserdova
#               miloslubov@gmail.com
#

set -euo pipefail
export LC_ALL=C

#Early check for build tools
cmake --version

cd $(dirname $(readlink -f "${BASH_SOURCE[0]}"))/..

mkdir -p ./build ./dist ./map_metrics
cp ./__init__.py ./map_metrics/__init__.py

cd ./build

NUMPROC=$(nproc)
echo "Running $NUMPROC parallel jobs"

LATEST=""

for PYBIN in /opt/python/cp3*/bin/python
do
    LATEST=${PYBIN}
    cmake -S .. -B . \
             -DPYTHON_EXECUTABLE:FILEPATH=${PYBIN} \
             -DCMAKE_BUILD_WITH_INSTALL_RPATH=TRUE \
             -DCMAKE_INSTALL_RPATH='$ORIGIN' \
             -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=$PWD/../bin \
             -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=$PWD/../map_metrics \
             -DCMAKE_CXX_FLAGS=-std=c++17 \
    && cmake --build .
done

cd ../
${LATEST} -m pip install $([[ -n "$VIRTUAL_ENV" ]] || echo "--user") -q build auditwheel
${LATEST} -m build --wheel --outdir ./dist/ .
auditwheel repair ./dist/*.whl
