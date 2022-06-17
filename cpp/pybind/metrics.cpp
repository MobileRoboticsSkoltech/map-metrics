// Copyright (c) 2022, Skolkovo Institute of Science and Technology (Skoltech)
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//
//  Created on: May 20, 2022
//       Author: Arthur Saliou
//               arthur.salio@gmail.com
//
#include <vector>
#include <Eigen/Core>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <metrics.h>
#include <orth_extract.h>
#include <config.h>

// Python: arg0: List[numpy.ndarray[numpy.float64[3, n]]], arg1: List[numpy.ndarray[numpy.float64[4, 4]]]) -> float

namespace py = pybind11;

std::vector<cilantro::VectorSet3d> extract_orthogonal(
       cilantro::VectorSet3d const & points , config::CustomConfig config, double eps=1e-1)
{
    return orth_extract::ExtractOrthogonalSubset(cilantro::PointCloud3d(points), config, eps);
}

PYBIND11_MODULE(map_metrics, m){
    m.doc() = R"(Map-Metrics library
A tool for evaluating the quality of odometry algorithm trajectory.)";

    py::module_ mcfg = m.def_submodule("config", "Config submodule");

    py::class_<config::CustomConfig>(mcfg, "CustomConfig")
            .def(py::init<int, double, int, int>());

    m.def("mpv", &metrics::GetMPV, "Mean Plane Variance trajectory metric",
          py::arg("pcs"),
          py::arg("poses"),
          py::arg("config") = config::CustomConfig());

    m.def("mme", &metrics::GetMME, "Mean Map Entropy trajectory metric",
          py::arg("pcs"),
          py::arg("poses"),
          py::arg("config") = config::CustomConfig());

    m.def("mom", &metrics::GetMOM, "Mutually Orthogonal Metric trajectory metric",
          py::arg("pcs"),
          py::arg("poses"),
          py::arg("config") = config::CustomConfig(),
          py::arg("orth_subset") = std::vector<cilantro::VectorSet3d>());

    m.def("extract_orthogonal", &extract_orthogonal, "Extract orthogonal plane subset from Point Cloud",
          py::arg("pc"),
          py::arg("config") = config::CustomConfig(),
          py::arg("eps") = 1e-1);
}