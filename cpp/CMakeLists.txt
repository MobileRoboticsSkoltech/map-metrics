set(SOURCES metrics.cpp metrics.h)

add_library(map_metrics SHARED ${SOURCES})

target_include_directories(map_metrics PUBLIC "${CMAKE_BINARY_DIR}/external/eigen/eigen/include/eigen3")
target_include_directories(map_metrics PUBLIC "${CMAKE_BINARY_DIR}/external/cilantro/cilantro/include")
target_include_directories(map_metrics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
