find_package(Boost 1.55.0 QUIET)

if (NOT Boost_FOUND)
    message("Couldn't find Boost >= 1.55 performing install")
    include(FetchContent)

    FetchContent_Declare(
        Boost
        GIT_REPOSITORY https://github.com/boostorg/boost.git
        GIT_TAG "boost-1.55.0"
    ) 

    FetchContent_MakeAvailable(Boost)
    include_directories(${Boost_INCLUDE_DIRS})
endif()
