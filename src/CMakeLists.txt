set(SOURCES
	Types.hpp
	Parameters.hpp

	auto_table.hpp
	auto_table.cpp

	BSPSubdivision.h
	BSPSubdivision.cpp

	Logger.hpp
	Logger.cpp

	Tetwild.h
	Tetwild.cpp
)

prepend_current_path(SOURCES)
tet_shell_copy_headers(${SOURCES})
tet_shell_set_source_group(${SOURCES})
target_sources(TetShell PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(external)
