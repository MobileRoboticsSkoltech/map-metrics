pybind11_add_module(map_metrics metrics.cpp)

find_package(Eigen3 3.4 REQUIRED)

target_link_libraries(map_metrics PRIVATE Eigen3::Eigen)
target_link_libraries(map_metrics PUBLIC map_metrics_cpp)
