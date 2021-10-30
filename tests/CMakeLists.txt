enable_testing()

add_subdirectory(gtest)

# Copying data from inputData to binary dir
file(GLOB PC_DATA ${CMAKE_SOURCE_DIR}/tests/data)
file(COPY ${PC_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(TEST_SOURCES main.cpp parse.cpp parse_utils.h)

add_executable(run_tests ${TEST_SOURCES})

# Including fetched content
target_include_directories(run_tests PUBLIC "${CMAKE_BINARY_DIR}/_deps/googletest-src/googletest/include")

target_link_libraries(run_tests PUBLIC gtest_main)
target_link_libraries(run_tests PUBLIC map_metrics)
target_link_libraries(run_tests PRIVATE Open3D::Open3D)

include(GoogleTest)
gtest_discover_tests(run_tests)