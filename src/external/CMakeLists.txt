set(SOURCES

	Rational.h
)

prepend_current_path(SOURCES)
tet_shell_copy_headers(${SOURCES})
tet_shell_set_source_group(${SOURCES})
target_sources(TetShell PRIVATE ${SOURCES})
