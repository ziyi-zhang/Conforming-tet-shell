# Distributed under the OSI-approved MIT License.  See accompanying
# file LICENSE or https://github.com/Crascit/DownloadProject for details.

cmake_minimum_required(VERSION 2.8.2)

project(cli11-download NONE)

include(ExternalProject)
ExternalProject_Add(cli11-download
                    URL;https://github.com/CLIUtils/CLI11/archive/v1.6.1.tar.gz;URL_MD5;48ef97262adb0b47a2f0a7edbda6e2aa
                    SOURCE_DIR          "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/cli11"
                    BINARY_DIR          "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/Release/cli11-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
)
