################################################################################
# Prepare dependencies
################################################################################

# Download external dependencies
include(TetShellDownloadExternal)

################################################################################
# Required libraries
################################################################################

# CL11
if(TET_SHELL_TOPLEVEL_PROJECT AND NOT TARGET CLI11::CLI11)
	tet_shell_download_cli11()
	add_subdirectory(${TET_SHELL_EXTERNAL}/cli11)
endif()

# fmt
if(NOT TARGET fmt::fmt)
	tet_shell_download_fmt()
	add_subdirectory(${TET_SHELL_EXTERNAL}/fmt)
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
	tet_shell_download_spdlog()

	# Create interface target
	add_library(spdlog INTERFACE)
	add_library(spdlog::spdlog ALIAS spdlog)
	target_include_directories(spdlog INTERFACE ${TET_SHELL_EXTERNAL}/spdlog/include)
	target_link_libraries(spdlog INTERFACE fmt::fmt)
	target_compile_definitions(spdlog INTERFACE SPDLOG_FMT_EXTERNAL)
endif()

# Libigl
if(NOT TARGET igl::core)
	tet_shell_download_libigl()

	# Import libigl targets
	list(APPEND CMAKE_MODULE_PATH "${TET_SHELL_EXTERNAL}/libigl/cmake")
	include(libigl)
endif()

# Geogram
if(NOT TARGET geogram::geogram)
	tet_shell_download_geogram()
	include(geogram)
endif()


# TBB
if(TET_SHELL_ENABLE_TBB AND NOT TARGET tbb::tbb)
	tet_shell_download_tbb()

	set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
	set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
	set(TBB_NO_DATE ON CACHE BOOL " " FORCE)

	add_subdirectory(${TET_SHELL_EXTERNAL}/tbb tbb)
	set_target_properties(tbb_static PROPERTIES
		INTERFACE_INCLUDE_DIRECTORIES "${TET_SHELL_EXTERNAL}/tbb/include"
	)
	if(NOT MSVC)
		set_target_properties(tbb_static PROPERTIES
			COMPILE_FLAGS "-Wno-implicit-fallthrough -Wno-missing-field-initializers -Wno-unused-parameter -Wno-keyword-macro"
		)
		set_target_properties(tbb_static PROPERTIES POSITION_INDEPENDENT_CODE ON)
	endif()
	add_library(tbb::tbb ALIAS tbb_static)
endif()

# C++11 threads
find_package(Threads REQUIRED)


# Json
if(NOT TARGET json)
	tet_shell_download_json()
	add_library(json INTERFACE)
	target_include_directories(json SYSTEM INTERFACE ${TET_SHELL_EXTERNAL}/json/include)
endif()



# winding number
tet_shell_download_windingnumber()
set(windingnumber_SOURCES
	${TET_SHELL_EXTERNAL}/windingnumber/SYS_Math.h
	${TET_SHELL_EXTERNAL}/windingnumber/SYS_Types.h
	${TET_SHELL_EXTERNAL}/windingnumber/UT_Array.cpp
	${TET_SHELL_EXTERNAL}/windingnumber/UT_Array.h
	${TET_SHELL_EXTERNAL}/windingnumber/UT_ArrayImpl.h
	${TET_SHELL_EXTERNAL}/windingnumber/UT_BVH.h
	${TET_SHELL_EXTERNAL}/windingnumber/UT_BVHImpl.h
	${TET_SHELL_EXTERNAL}/windingnumber/UT_FixedVector.h
	${TET_SHELL_EXTERNAL}/windingnumber/UT_ParallelUtil.h
	${TET_SHELL_EXTERNAL}/windingnumber/UT_SmallArray.h
	${TET_SHELL_EXTERNAL}/windingnumber/UT_SolidAngle.cpp
	${TET_SHELL_EXTERNAL}/windingnumber/UT_SolidAngle.h
	${TET_SHELL_EXTERNAL}/windingnumber/VM_SIMD.h
	${TET_SHELL_EXTERNAL}/windingnumber/VM_SSEFunc.h
)

add_library(fast_winding_number ${windingnumber_SOURCES})
target_link_libraries(fast_winding_number PRIVATE tbb::tbb)
target_compile_features(fast_winding_number PRIVATE ${CXX17_FEATURES})
target_include_directories(fast_winding_number PUBLIC "${TET_SHELL_EXTERNAL}/")
