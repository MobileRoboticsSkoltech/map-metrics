cmake_minimum_required(VERSION 3.11)

project(map-metrics)

# Compiler prerequisites
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # # # # # OPTIONS # # # # # #
option(GLIBCXX_USE_CXX11_ABI   "Set -D_GLIBCXX_USE_CXX11_ABI=1"       OFF)
# Build Tests
option(BUILD_TESTS "Build an executable with tests" OFF)
# Python Bindings
option(BUILD_PYTHON_MODULE "Build the Python library interface" ON)
# Address sanitizer
option(USE_ASAN "Address sanitizer" OFF) 
# Memory sanitizer
option(USE_MSAN "Memory sanitizer" OFF)

if (WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    set(BUILD_SHARED_LIBS TRUE)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

# By default select Release build
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No CMAKE_BUILD_TYPE specified, default to Release.")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# External dependencies
add_subdirectory(external)
# Map-Metrics sources
add_subdirectory(cpp)
