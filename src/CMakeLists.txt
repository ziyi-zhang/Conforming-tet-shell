set(SOURCES
	Types.hpp
	Parameters.hpp

	auto_table.hpp
	auto_table.cpp

	Logger.hpp
	Logger.cpp
)

prepend_current_path(SOURCES)
tet_shell_copy_headers(${SOURCES})
tet_shell_set_source_group(${SOURCES})
target_sources(TetShell PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(external)
