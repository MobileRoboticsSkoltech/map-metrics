set(SOURCES metrics.cpp metrics.h metrics_utils.cpp metrics_utils.h)

add_library(map_metrics SHARED ${SOURCES})

target_include_directories(map_metrics PUBLIC "${CMAKE_BINARY_DIR}/_deps/eigen-src")
target_include_directories(map_metrics PUBLIC "${CMAKE_BINARY_DIR}/external/cilantro/cilantro/include")
target_include_directories(map_metrics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
