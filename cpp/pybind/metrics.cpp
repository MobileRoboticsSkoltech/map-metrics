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
double py_mpv(std::vector<cilantro::VectorSet3d> const & pcs_points, std::vector<Eigen::Matrix4d> const & ts){
    std::vector<cilantro::PointCloud3d> pcs(pcs_points.size());
    for (size_t i = 0; i < pcs_points.size(); ++i){
        pcs[i] = cilantro::PointCloud3d(pcs_points[i]);
    }

    return map_metrics::mpv(pcs, ts);
}

double py_mme(std::vector<cilantro::VectorSet3d> const & pcs_points, std::vector<Eigen::Matrix4d> const & ts){
    std::vector<cilantro::PointCloud3d> pcs(pcs_points.size());
    for (size_t i = 0; i < pcs_points.size(); ++i){
        pcs[i] = cilantro::PointCloud3d(pcs_points[i]);
    }

    return map_metrics::mme(pcs, ts);
}

namespace py = pybind11;

PYBIND11_MODULE(map_metrics_py, m){
    m.doc() = "Baseline of MPV and MME metrics";
    m.def("mpv", py::overload_cast<
          const std::vector<cilantro::VectorSet3d> &,
          const std::vector<Eigen::Matrix4d> &>(&py_mpv));
    m.def("mme", py::overload_cast<
          const std::vector<cilantro::VectorSet3d> &,
          const std::vector<Eigen::Matrix4d> &>(&py_mme));
}
