include(FetchContent)

FetchContent_Declare(
        cilantro
        GIT_REPOSITORY https://github.com/kzampog/cilantro.git
        GIT_TAG origin/master
)

FetchContent_GetProperties(cilantro)
if(NOT cilantro_POPULATED)
    FetchContent_Populate(cilantro)
    set(cilantro_BUILD_EXAMPLES OFF)
    target_include_directories(map_metrics PUBLIC "${CMAKE_BINARY_DIR}/_deps/cilantro-src/include")
endif()

# TODO: Cilantro root CMakeLists.txt:70 somehow doesn't add config.hpp file, we do it manually
configure_file("${CMAKE_BINARY_DIR}/_deps/cilantro-src/config.hpp.in"
        "${CMAKE_BINARY_DIR}/_deps/cilantro-src/include/cilantro/config.hpp")
