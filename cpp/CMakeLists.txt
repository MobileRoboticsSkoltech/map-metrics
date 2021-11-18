set(SOURCES metrics.cpp metrics.h)

add_library(map_metrics SHARED ${SOURCES})

target_include_directories(map_metrics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
