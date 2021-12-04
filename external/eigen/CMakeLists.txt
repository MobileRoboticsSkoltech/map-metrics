include(ExternalProject)
ExternalProject_Add(
        external_eigen
        PREFIX eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4
        CMAKE_ARGS
        -DEIGEN_BUILD_DOC=OFF
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)
