find_package(eigen QUIET)

if (NOT eigen_FOUND)
    include(FetchContent)

    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4
    )
    set(eigen_EIGEN_BUILD_DOC OFF)
    set(eigen_BUILD_TESTS OFF)

    FetchContent_MakeAvailable(eigen)
endif()
