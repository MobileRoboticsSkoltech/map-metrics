include(FetchContent)

FetchContent_Declare(
        eigen 
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4
)

# FetchContent_GetProperties(eigen)
# if (NOT eigen_POPULATED)
#     FetchContent_Populate(eigen)
#     set(EIGEN_BUILD_DOC OFF)
#     set(EIGEN_BUILD_TESTS OFF)
#     add_subdirectory(${eigen_SOURCE_DIR} ${eigen_BINARY_DIR})
# endif()

FetchContent_MakeAvailable(eigen)
