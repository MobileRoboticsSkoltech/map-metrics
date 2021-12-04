pybind11_add_module(map_metrics_py metrics.cpp)

target_link_libraries(map_metrics_py PUBLIC map_metrics)
