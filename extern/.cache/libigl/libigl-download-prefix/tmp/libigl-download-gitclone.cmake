
if(NOT "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/.cache/libigl/libigl-download-prefix/src/libigl-download-stamp/libigl-download-gitinfo.txt" IS_NEWER_THAN "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/.cache/libigl/libigl-download-prefix/src/libigl-download-stamp/libigl-download-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/.cache/libigl/libigl-download-prefix/src/libigl-download-stamp/libigl-download-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/libigl"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/libigl'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"  clone  "https://github.com/libigl/libigl.git" "libigl"
    WORKING_DIRECTORY "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/libigl/libigl.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"  checkout v2.2.0 --
  WORKING_DIRECTORY "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/libigl"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'v2.2.0'")
endif()

execute_process(
  COMMAND "/usr/bin/git"  submodule update --recursive --init 
  WORKING_DIRECTORY "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/libigl"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/libigl'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/.cache/libigl/libigl-download-prefix/src/libigl-download-stamp/libigl-download-gitinfo.txt"
    "/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/.cache/libigl/libigl-download-prefix/src/libigl-download-stamp/libigl-download-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/Users/ziyizhang/Projects/ShellMesh/Conforming-tet-shell/extern/.cache/libigl/libigl-download-prefix/src/libigl-download-stamp/libigl-download-gitclone-lastrun.txt'")
endif()

