pybind11_add_module(map_metrics metrics.cpp)

find_package(Eigen3 3.4 REQUIRED)

# TODO: remove?
target_link_libraries(map_metrics PRIVATE ${cilantro_LIBRARIES})
target_include_directories(map_metrics PRIVATE "${CMAKE_BINARY_DIR}/external/cilantro/cilantro/include")

target_link_libraries(map_metrics PRIVATE Eigen3::Eigen)
target_link_libraries(map_metrics PUBLIC map_metrics_cpp)
