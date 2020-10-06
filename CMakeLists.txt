################################################################################
cmake_minimum_required(VERSION 3.8)
project(TetShell)
################################################################################

# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
	set(TET_SHELL_TOPLEVEL_PROJECT OFF)
else()
	set(TET_SHELL_TOPLEVEL_PROJECT ON)
endif()

if(INPUT_THIRD_PARTY_DIR)
	set(TET_SHELL_EXTERNAL ${CMAKE_CURRENT_SOURCE_DIR}/${INPUT_THIRD_PARTY_DIR}/)
else()
	set(TET_SHELL_EXTERNAL ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Color output
include(UseColors)

# Prepend function
include(PrependCurrentPath)

# TetShell utils
include(TetShellUtils)

# Extra warnings
include(Warnings)

# Sort projects inside the solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################

# if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
# 	foreach(config ${CMAKE_CONFIGURATION_TYPES})
# 		string(TOUPPER ${config} config)
# 		string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
# 		string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
# 	endforeach()
# endif()

################################################################################

# TetShell options
option(TET_SHELL_ENABLE_TBB         "Enable TBB"                    ON)
option(TET_SHELL_USE_FLOAT          "Use floats instead of double"  OFF)

# Options for libigl modules
option(LIBIGL_USE_STATIC_LIBRARY          "Use libigl as static library" OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
option(LIBIGL_WITH_PNG               "Use PNG"                      OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         OFF)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)

# Setup dependencies
include(TetShellDependencies)

################################################################################
# TetShell library
################################################################################

find_package(GMPfTetWild)
IF(NOT ${GMP_FOUND})
	MESSAGE(FATAL_ERROR "Cannot find GMP")
ENDIF()

# find_package(MPFR)
# IF(NOT ${MPFR_FOUND})
# 	MESSAGE(FATAL_ERROR "Cannot find MPFR")
# ENDIF()

# add_library() can only be called without any source since CMake 3.11 ...
add_library(${PROJECT_NAME} src/Logger.cpp)

# Public include directory for TetShell
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR}/include)

# set(MESH_PATH "${CMAKE_CURRENT_SOURCE_DIR}/tests/")
# target_compile_definitions(${PROJECT_NAME} PUBLIC -DTET_SHELL_MESH_PATH=\"${MESH_PATH}\")

# Extra warnings
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)

# Use C++11
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
# if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
# 	target_compile_options(${PROJECT_NAME} PUBLIC "/Zc:__cplusplus")
# endif()
# target_compile_definitions(${PROJECT_NAME} PUBLIC -DEIGEN_STACK_ALLOCATION_LIMIT=8388608)

################################################################################
# Required libraries
################################################################################

if(TET_SHELL_USE_FLOAT)
	target_compile_definitions(${PROJECT_NAME} PUBLIC -DTET_SHELL_USE_FLOAT)
endif()

target_include_directories(${PROJECT_NAME} PUBLIC ${GMP_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME}
	PUBLIC
		igl::core
		geogram::geogram
		spdlog::spdlog
		Threads::Threads
		fast_winding_number
		json
		${GMP_LIBRARIES}
		# ${MPFR_LIBRARIES}
)
if(TET_SHELL_ENABLE_TBB)
	target_link_libraries(${PROJECT_NAME} PUBLIC tbb::tbb)
	target_compile_definitions(${PROJECT_NAME} PUBLIC TET_SHELL_USE_TBB)
endif()

################################################################################
# TetShell binary
################################################################################

# Main executable
if(TET_SHELL_TOPLEVEL_PROJECT)
	add_executable(${PROJECT_NAME}_bin src/main.cpp)
	target_compile_features(${PROJECT_NAME}_bin PUBLIC ${CXX14_FEATURES})

	target_link_libraries(${PROJECT_NAME}_bin
		PUBLIC
			${PROJECT_NAME}
			CLI11::CLI11
			warnings::all
	)

	if(TARGET igl::tetgen)
		target_link_libraries(${PROJECT_NAME} PUBLIC igl::tetgen)
		target_compile_definitions(${PROJECT_NAME} PUBLIC -DLIBIGL_WITH_TETGEN)
	endif()

	if(NOT (${CMAKE_VERSION} VERSION_LESS "3.6.0"))
		set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME}_bin)
	endif()
endif()

################################################################################
# Subdirectories
################################################################################

# Other sources
add_subdirectory(src)

# Compile extras only if this is a top-level project
if(TET_SHELL_TOPLEVEL_PROJECT)
	# Unit tests
	include(CTest)
	enable_testing()
	add_subdirectory(tests)
endif()
