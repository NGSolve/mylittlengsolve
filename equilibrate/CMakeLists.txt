add_ngsolve_python_module(libequilibrate equilibrate.cpp)


set_target_properties(libequilibrate PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/${NETGEN_PYTHON_RPATH}")
install(TARGETS libequilibrate DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON})
