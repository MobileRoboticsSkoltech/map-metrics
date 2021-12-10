include(ExternalProject)
ExternalProject_Add(
    external_cilantro
    PREFIX cilantro
    GIT_REPOSITORY https://github.com/kzampog/cilantro.git
    GIT_TAG origin/master
    CMAKE_ARGS
        -DBUILD_EXAMPLES=OFF
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DBUILD_SHARED_LIBS=OFF
)
