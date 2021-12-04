include(ExternalProject)
ExternalProject_Add(
        external_eigen
        PREFIX eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4
        CMAKE_ARGS
        -DEIGEN_BUILD_DOC=OFF
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)

#include(FetchContent)
#
#FetchContent_Declare(
#        eigen
#        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
#        GIT_TAG 3.4
#)
#
#FetchContent_GetProperties(eigen)
#if(NOT eigen_POPULATED)
#    FetchContent_Populate(eigen)
#    set(eigen_EIGEN_BUILD_DOC OFF)
#    set(eigen_BUILD_TESTS OFF)
#    add_subdirectory(${eigen_SOURCE_DIR} ${eigen_BINARY_DIR})
#    # target_include_directories(map_metrics PUBLIC "${CMAKE_BINARY_DIR}/_deps/eigen-src")
#endif()
