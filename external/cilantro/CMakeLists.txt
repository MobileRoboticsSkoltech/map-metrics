find_package(cilantro QUIET)

if (NOT cilantro_FOUND)
    include(ExternalProject)
    ExternalProject_Add(
        cilantro_external
        PREFIX cilantro
        GIT_REPOSITORY https://github.com/kzampog/cilantro.git
        GIT_TAG origin/master
        CMAKE_ARGS
            -DBUILD_EXAMPLES=OFF
            -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
            # -DBUILD_SHARED_LIBS=OFF
    )
    ExternalProject_Get_Property(cilantro_external INSTALL_DIR)
    set(cilantro_LIBRARIES 
    "${INSTALL_DIR}/lib/${CMAKE_SHARED_MODULE_PREFIX}cilantro${CMAKE_SHARED_LIBRARY_SUFFIX}" 
    CACHE INTERNAL "" FORCE
    )
endif()
