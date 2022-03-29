cmake_minimum_required(VERSION 3.19.2)
project(map-metrics)

# Compiler prerequisites
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # # # # # OPTIONS # # # # # #
option(GLIBCXX_USE_CXX11_ABI   "Set -D_GLIBCXX_USE_CXX11_ABI=1"       OFF)
# Set this flag ON if you want to build tests
option(BUILD_TESTS "Build an executable with tests" ON)
# Set this flag ON if you want to build Python library interface
option(BUILD_PYTHON_MODULE "Build the python module" OFF)
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
    # TODO: Set Release to default 
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# # Enable sanitizer by default in Debug mode
# if (CMAKE_BUILD_TYPE MATCHES Debug)
#     set(USE_ASAN ON)
# endif()

# If Build type is Release, then turn on optimizations
if (CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
endif()

if (USE_ASAN)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
endif()

if (USE_MSAN)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory -fno-omit-frame-pointer")
endif()

# Prevent In-Source Builds
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
message(FATAL_ERROR "In-source build detected!")
endif()

# External dependencies
add_subdirectory(external)
# Be your own benchmark library
add_subdirectory(cpp)
