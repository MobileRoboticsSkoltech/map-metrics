pybind11_add_module(map_metrics metrics.cpp)

target_link_libraries(map_metrics PUBLIC map_metrics_cpp)
