cmake_minimum_required(VERSION 3.3)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

################################################################################


################################################################################
# Options
################################################################################
# tetshell
option(TETSHELL_LIBONLY "Use TetShell as a library, without executables" OFF)
# tetwild
option(TETWILD_WITH_HUNTER "Use Hunter to download and configure Boost" OFF)
option(TETWILD_WITH_ISPC   "Use ISPC"                                   OFF)
option(TETWILD_WITH_ASAN   "Use ASAN" 									OFF)
# libigl library
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_CGAL             "Use CGAL"           ON)
if (NOT TETSHELL_LIBONLY)
option(LIBIGL_WITH_TETGEN           "Use Tetgen"         ON)
endif()
# gui
option(LIBIGL_WITH_EMBREE           "Use Embree"         OFF)
option(LIBIGL_WITH_NANOGUI          "Use Nanogui menu"   OFF)
option(LIBIGL_WITH_OPENGL           "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW      "Use GLFW"           OFF)
option(LIBIGL_WITH_VIEWER           "Use OpenGL viewer"  OFF)

if(TETWILD_WITH_HUNTER)
	# Needs to be set before the main project... argh =/
	include(HunterGate)
	HunterGate(
		URL "https://github.com/ruslo/hunter/archive/v0.23.25.tar.gz"
		SHA1 "cb75cce9a3a8d552e70e7118f3203eb4ac05c201"
	)
endif()

################################################################################
# Project name
################################################################################

project(ConformingTetShell)

################################################################################
# Settings
################################################################################
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(TETWILD_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/extern")

# Color output
include(UseColors)

# Use folder in Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Extra warnings
include(Warnings)

# Export compile flags (used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Generate position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# 3rd party libraries
################################################################################
include(TetWildDependencies)

################################################################################
# TetWild
################################################################################
# Build static library for executable
add_library(libTetShell 
		include/tetwild/Args.h
		include/tetwild/Exception.h
		include/tetwild/Logger.h
		include/tetwild/tetwild.h
		src/tetwild/BSPSubdivision.cpp
		src/tetwild/BSPSubdivision.h
		src/tetwild/CGALTypes.h
		src/tetwild/Common.cpp
		src/tetwild/Common.h
		src/tetwild/DelaunayTetrahedralization.cpp
		src/tetwild/DelaunayTetrahedralization.h
		src/tetwild/DistanceQuery.cpp
		src/tetwild/DistanceQuery.h
		src/tetwild/EdgeCollapser.cpp
		src/tetwild/EdgeCollapser.h
		src/tetwild/EdgeRemover.cpp
		src/tetwild/EdgeRemover.h
		src/tetwild/EdgeSplitter.cpp
		src/tetwild/EdgeSplitter.h
		src/tetwild/ForwardDecls.h
		src/tetwild/InoutFiltering.cpp
		src/tetwild/InoutFiltering.h
		src/tetwild/LocalOperations.cpp
		src/tetwild/LocalOperations.h
		src/tetwild/Logger.cpp
		src/tetwild/MeshConformer.cpp
		src/tetwild/MeshConformer.h
		src/tetwild/MeshRefinement.cpp
		src/tetwild/MeshRefinement.h
		src/tetwild/Preprocess.cpp
		src/tetwild/Preprocess.h
		src/tetwild/SimpleTetrahedralization.cpp
		src/tetwild/SimpleTetrahedralization.h
		src/tetwild/State.cpp
		src/tetwild/State.h
		src/tetwild/TetmeshElements.cpp
		src/tetwild/TetmeshElements.h
		src/tetwild/tetwild.cpp
		src/tetwild/VertexSmoother.cpp
		src/tetwild/VertexSmoother.h
		src/tetwild/geogram/mesh_AABB.cpp
		src/tetwild/geogram/mesh_AABB.h

		src/shell/common.hpp
		src/shell/Label.cpp
		src/shell/Label.h
		src/shell/Shell.cpp
		src/shell/Shell.h
		src/shell/Utils.cpp
		src/shell/Utils.h
		src/shell/ShellCheck.cpp
		src/shell/ShellCheck.h
		src/shell/TetMeshCheck.cpp
		src/shell/TetMeshCheck.h
		src/shell/OptimalEnergy.cpp
		src/shell/OptimalEnergy.h
)
target_include_directories(libTetShell
	PUBLIC
		src
		include
)
target_link_libraries(libTetShell
	PUBLIC
		geogram
		igl::core
		pymesh::pymesh
		spdlog::spdlog
		fmt::fmt
		highfive
	PRIVATE
		igl::cgal
		warnings::all
)
set_target_properties(libTetShell PROPERTIES OUTPUT_NAME "tetwild")

# ispc
if (TETWILD_WITH_ISPC)
	message(STATUS "Compiling energy with ISPC")
	add_subdirectory(src/ispc)
	ispc_add_energy(libTetShell)
endif()

# ASAN
if (TETWILD_WITH_ASAN)
	message(STATUS "Compiling with ASAN")
	target_compile_options(libTetShell PUBLIC "-fsanitize=address")
	target_link_options(libTetShell PUBLIC "-fsanitize=address")
endif()

# Building executable test
if (NOT TETSHELL_LIBONLY)
add_executable(TetShellTest EXCLUDE_FROM_ALL src/test/main_test.cpp)
target_link_libraries(TetShellTest
		libTetShell
		CLI11::CLI11
		igl::cgal
		warnings::all
)
target_include_directories(TetShellTest PRIVATE src)
igl_copy_cgal_dll(TetShellTest)

# Building executable tetgen
if (TETSHELL_TEST_TETGEN)
	add_executable(TetGen EXCLUDE_FROM_ALL src/test/main_tetgen.cpp)
	target_link_libraries(TetGen
			libTetShell
			igl::core
			igl::tetgen
			warnings::all
	)
	target_include_directories(TetGen PRIVATE src)
endif()

# Building main executable
add_executable(TetShell src/main.cpp)
target_link_libraries(TetShell
		libTetShell
		CLI11::CLI11
		igl::cgal
		warnings::all
)
target_include_directories(TetShell PRIVATE src)
igl_copy_cgal_dll(TetShell)

# Install
install(TARGETS TetShell RUNTIME DESTINATION bin)
endif()

################################################################################
# Folders for Visual Studio/XCode IDEs
################################################################################

# geogram
set_target_properties(geogram PROPERTIES FOLDER extern/geogram)
set_target_properties(geogram_third_party PROPERTIES FOLDER extern/geogram)
set_target_properties(uninstall PROPERTIES FOLDER extern/geogram)
# pymesh
set_target_properties(pymesh_tiny PROPERTIES FOLDER extern/pymesh)
