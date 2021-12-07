pybind11_add_module(map_metrics_cpp metrics.cpp)

target_link_libraries(map_metrics_cpp PUBLIC map_metrics)
