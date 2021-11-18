cmake_minimum_required(VERSION 3.19.2)
project(map-metrics)

# Compiler prerequisites
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# # # # # # OPTIONS # # # # # #
option(GLIBCXX_USE_CXX11_ABI   "Set -D_GLIBCXX_USE_CXX11_ABI=1"       OFF)
option(STATIC_WINDOWS_RUNTIME  "Use static (MT/MTd) Windows runtime"  ON )
# Set this flag ON if you want to build tests
option(BUILD_TESTS "Build an executable with tests" OFF)

# By default select Debug build
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No CMAKE_BUILD_TYPE specified, default to Debug.")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# If Build type is Release, then turn on optimizations
if (CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
endif()

# Prevent In-Source Builds
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
message(FATAL_ERROR "In-source build detected!")
endif()

# Be your own benchmark library
add_subdirectory(cpp)
# External dependencies
add_subdirectory(external)

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()
