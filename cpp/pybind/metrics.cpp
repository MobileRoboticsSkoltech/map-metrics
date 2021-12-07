//
// Created by achains on 04.12.2021.
//
#include <vector>
#include <Eigen/Core>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <metrics.h>

// Python: arg0: List[numpy.ndarray[numpy.float64[3, n]]], arg1: List[numpy.ndarray[numpy.float64[4, 4]]]) -> float

namespace py = pybind11;

PYBIND11_MODULE(map_metrics_py, m){
    m.doc() = "Baseline of MPV and MME metrics";
    m.def("mpv", py::overload_cast<
          const std::vector<cilantro::VectorSet3d> &,
          const std::vector<Eigen::Matrix4d> &, int, double>(&metrics::GetMPV));
    m.def("mme", py::overload_cast<
          const std::vector<cilantro::VectorSet3d> &,
          const std::vector<Eigen::Matrix4d> &, int, double>(&metrics::GetMME));
}