# Distributed under the OSI-approved MIT License.  See accompanying
# file LICENSE or https://github.com/Crascit/DownloadProject for details.

cmake_minimum_required(VERSION 2.8.2)

project(libigl-download NONE)

include(ExternalProject)
ExternalProject_Add(libigl-download
                    GIT_REPOSITORY;https://github.com/libigl/libigl.git;GIT_TAG;v2.2.0
                    SOURCE_DIR          "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/libigl"
                    BINARY_DIR          "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/Release/libigl-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
)
