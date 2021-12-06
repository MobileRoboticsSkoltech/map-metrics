include(ExternalProject)
ExternalProject_Add(
    external_cilantro
    PREFIX cilantro
    GIT_REPOSITORY https://github.com/kzampog/cilantro.git
    GIT_TAG origin/master
    CMAKE_ARGS
        -DBUILD_EXAMPLES=OFF
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)
