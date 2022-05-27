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
#  build-wheels-macOS.sh
# 
#  Created on: Mar 22, 2021
#       Author: Lyubov Miloserdova
#               miloslubov@gmail.com
#

set -euo pipefail
export LC_ALL=C
export MACOSX_DEPLOYMENT_TARGET=10.15

cd $(dirname $(greadlink -f "${BASH_SOURCE[0]}"))/..
mkdir -p ./build ./dist ./map_metrics
cp ./__init__.py ./map_metrics/__init__.py

cd ./build

export DYLD_LIBRARY_PATH=${BOOST_ROOT}/lib

for PYBIN in /Users/runner/hostedtoolcache/Python/3.*/x64/bin/python*?[0-9]
do
    cmake .. -DPYTHON_EXECUTABLE:FILEPATH=$PYBIN \
             -DCMAKE_BUILD_WITH_INSTALL_RPATH=TRUE \
             -DCMAKE_INSTALL_RPATH="@loader_path" \
             -DBOOST_ROOT=${BOOST_ROOT} \
             -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=$PWD/../bin \
             -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=$PWD/../map_metrics \
             -DCMAKE_CXX_FLAGS="-std=c++1z" \
    && cmake --build .
done

cd ../
python3 -m pip install --user -q build
python3 -m build --wheel --outdir dist/ .
