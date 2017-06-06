add_ngsolve_python_module(myngspy myngspy.cpp)

set_target_properties(myngspy PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/${NETGEN_PYTHON_RPATH}")
target_link_libraries(myngspy myngsolve)
install(TARGETS myngspy DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON})
