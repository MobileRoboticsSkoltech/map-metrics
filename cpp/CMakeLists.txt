set(SOURCES metrics.cpp metrics.h metrics_utils.cpp orth_extract.cpp 
metrics_utils.h config.h orth_extract.h clustering.h clustering.cpp
max_clique.h max_clique.cpp)

add_library(map_metrics_cpp SHARED ${SOURCES})

add_dependencies(map_metrics_cpp external_cilantro)

target_include_directories(map_metrics_cpp PUBLIC "${CMAKE_BINARY_DIR}/_deps/eigen-src")
target_include_directories(map_metrics_cpp PUBLIC "${CMAKE_BINARY_DIR}/external/cilantro/cilantro/include")
target_include_directories(map_metrics_cpp PUBLIC "${CMAKE_BINARY_DIR}/external/alglib/include")
target_include_directories(map_metrics_cpp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(map_metrics_cpp alglib)

if (BUILD_PYTHON_MODULE)
    add_subdirectory("${CMAKE_SOURCE_DIR}/external/pybind11" "${CMAKE_BINARY_DIR}/external/pybind11")
    add_subdirectory(pybind)
endif()

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()