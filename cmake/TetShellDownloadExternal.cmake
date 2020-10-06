################################################################################

include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(TET_SHELL_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(TET_SHELL_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(tet_shell_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${TET_SHELL_EXTERNAL}/${name}
        DOWNLOAD_DIR ${TET_SHELL_EXTERNAL}/.cache/${name}
        QUIET
        ${TET_SHELL_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

## libigl
function(tet_shell_download_libigl)
    tet_shell_download_project(libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG        015ac35cd135e3a40c3e71443eeda4a9a4ccf7f5
    )
endfunction()

## Json
function(tet_shell_download_json)
    tet_shell_download_project(json
        GIT_REPOSITORY https://github.com/jdumas/json
        GIT_TAG        0901d33bf6e7dfe6f70fd9d142c8f5c6695c6c5b
    )
endfunction()

## Catch2
function(tet_shell_download_catch2)
    tet_shell_download_project(Catch2
        URL          https://github.com/catchorg/Catch2/archive/v2.4.2.tar.gz
        URL_MD5      26927b878b1f42633f15a9ef1c4bd8e7
    )
endfunction()

## CLI11
function(tet_shell_download_cli11)
    tet_shell_download_project(cli11
            URL     https://github.com/CLIUtils/CLI11/archive/v1.8.0.tar.gz
            URL_MD5 5e5470abcb76422360409297bfc446ac
    )
endfunction()

## tbb
function(tet_shell_download_tbb)
    tet_shell_download_project(tbb
        GIT_REPOSITORY https://github.com/wjakob/tbb.git
        GIT_TAG        ddbe45cd3ad89df9a84cd77013d5898fc48b8e89
    )
endfunction()

## fmt
function(tet_shell_download_fmt)
    tet_shell_download_project(fmt
        URL      https://github.com/fmtlib/fmt/archive/5.3.0.tar.gz
        URL_MD5  1015bf3ff2a140dfe03de50ee2469401
    )
endfunction()

## spdlog
function(tet_shell_download_spdlog)
    tet_shell_download_project(spdlog
        URL         https://github.com/gabime/spdlog/archive/v1.3.1.tar.gz
        URL_MD5     3c17dd6983de2a66eca8b5a0b213d29f
    )
endfunction()

## Geogram LGPL
function(tet_shell_download_geogram)
    tet_shell_download_project(geogram
        GIT_REPOSITORY   https://github.com/polyfem/geogram.git
        GIT_TAG          e6b9612f1146370e40deaa341b4dd7ef90502102
    )
endfunction()

## aabbcc
function(tet_shell_download_aabbcc)
    tet_shell_download_project(aabbcc
            GIT_REPOSITORY https://github.com/lohedges/aabbcc.git
            GIT_TAG        0c85e61362d384d70c71946826bfed0fb24a74ba
            )
endfunction()

## winding number
function(tet_shell_download_windingnumber)
    tet_shell_download_project(windingnumber
            GIT_REPOSITORY https://github.com/alecjacobson/WindingNumber.git
            GIT_TAG        1e6081e52905575d8e98fb8b7c0921274a18752f
            )
endfunction()
