set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TEST_DIR}/tcl )
add_library( tcldemo SHARED ../../tcl/tcldemo.cpp )
target_link_libraries( tcldemo ngsolve )
test_pde( tcldemo.pde ${TEST_DIR}/tcl ${TEST_DIR}/tcl/libtcldemo.so )
