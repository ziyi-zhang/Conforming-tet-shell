# Distributed under the OSI-approved MIT License.  See accompanying
# file LICENSE or https://github.com/Crascit/DownloadProject for details.

cmake_minimum_required(VERSION 2.8.2)

project(fmt-download NONE)

include(ExternalProject)
ExternalProject_Add(fmt-download
                    GIT_REPOSITORY;https://github.com/fmtlib/fmt.git;GIT_TAG;5.2.0
                    SOURCE_DIR          "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/fmt"
                    BINARY_DIR          "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/Release/fmt-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
)
