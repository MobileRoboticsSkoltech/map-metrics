enable_testing()

add_subdirectory(gtest)

# Copying data from inputData to binary dir
file(GLOB PC_DATA ${CMAKE_SOURCE_DIR}/cpp/tests/data)
file(COPY ${PC_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(run_tests 
    main.cpp
    test_parse.cpp 
    test_base_metrics.cpp
    test_depth.cpp
    test_orthogonal.cpp
    parse_utils.h
    parse_utils.cpp
    pc_utils.h
    pc_utils.cpp
)

find_package(Eigen3 3.4 REQUIRED)
if (Eigen3_FOUND)
    target_link_libraries(run_tests PRIVATE Eigen3::Eigen)
endif()

target_link_libraries(run_tests PRIVATE ${cilantro_LIBRARIES})
target_include_directories(run_tests PRIVATE "${CMAKE_BINARY_DIR}/external/cilantro/cilantro/include")


target_link_libraries(run_tests PRIVATE gtest_main)
target_link_libraries(run_tests PRIVATE map_metrics_cpp)

target_compile_definitions(run_tests PRIVATE 
    TEST_ORTHOGONAL_METRICS 
    TEST_BASE_DEPTH_METRICS
)

include(GoogleTest)
gtest_discover_tests(run_tests)