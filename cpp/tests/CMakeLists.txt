enable_testing()

add_subdirectory(gtest)

# Copying data from inputData to binary dir
file(GLOB PC_DATA ${CMAKE_SOURCE_DIR}/cpp/tests/data)
file(COPY ${PC_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(TEST_SOURCES main.cpp test_parse.cpp test_base_metrics.cpp parse_utils.cpp parse_utils.h pc_utils.h)

add_executable(run_tests ${TEST_SOURCES})

# Including fetched content
target_include_directories(run_tests PUBLIC "${CMAKE_BINARY_DIR}/_deps/googletest-src/googletest/include")

target_link_libraries(run_tests PUBLIC gtest_main)
target_link_libraries(run_tests PUBLIC map_metrics_cpp)

include(GoogleTest)
gtest_discover_tests(run_tests)