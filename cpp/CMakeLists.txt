add_library(map_metrics_cpp SHARED 
  metrics.h
  metrics.cpp
  metrics_utils.h
  metrics_utils.cpp 
  orth_extract.h
  orth_extract.cpp
  clustering.h 
  clustering.cpp
  max_clique.h 
  max_clique.cpp
)

if (CMAKE_BUILD_TYPE MATCHES Release)
  target_compile_options(map_metrics_cpp PRIVATE "-O2")
endif()

find_package(Boost 1.55.0 REQUIRED COMPONENTS graph)
if(Boost_FOUND)
  target_link_libraries(map_metrics_cpp PUBLIC Boost::graph)
  # target_include_directories(map_metrics_cpp PRIVATE "${Boost_INCLUDE_DIRS}")
endif()

find_package(Eigen3 3.4 REQUIRED)
if (Eigen3_FOUND)
  target_link_libraries(map_metrics_cpp PRIVATE Eigen3::Eigen)
endif()

target_include_directories(map_metrics_cpp PRIVATE "${CMAKE_BINARY_DIR}/external/cilantro/cilantro/include")
# TODO: Is it really necessary to specify target_include_dir for alglib here?
target_include_directories(map_metrics_cpp PUBLIC "${CMAKE_BINARY_DIR}/external/alglib/include")
target_include_directories(map_metrics_cpp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(map_metrics_cpp PRIVATE alglib)
target_link_libraries(map_metrics_cpp PRIVATE ${cilantro_LIBRARIES})

if (BUILD_PYTHON_MODULE)
    add_subdirectory("${CMAKE_SOURCE_DIR}/external/pybind11" "${CMAKE_BINARY_DIR}/external/pybind11")
    add_subdirectory(pybind)
endif()

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()