//
// Created by achains on 04.12.2021.
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
    m.doc() = "Baseline of MPV and MME metrics";
    py::module_ mcfg = m.def_submodule("config", "Config submodule");

    py::class_<config::CustomConfig>(mcfg, "CustomConfig")
            .def(py::init<int, double, int, int>());

    m.def("mpv", py::overload_cast<
            std::vector<cilantro::VectorSet3d> const &,
            std::vector<Eigen::Matrix4d> const &,
            config::CustomConfig
                >(&metrics::GetMPV));

    m.def("mme", py::overload_cast<
            std::vector<cilantro::VectorSet3d> const &,
            std::vector<Eigen::Matrix4d> const &,
            config::CustomConfig
                >(&metrics::GetMME));

    m.def("mom", py::overload_cast<
            std::vector<cilantro::VectorSet3d> const &,
            std::vector<Eigen::Matrix4d> const &,
            config::CustomConfig,
            std::vector<cilantro::VectorSet3d> const &
                >(&metrics::GetMOM));

    m.def("extract_orthogonal", py::overload_cast<
            cilantro::VectorSet3d const &,
            config::CustomConfig,
            double>(&extract_orthogonal));
}