set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TEST_DIR}/precond )
add_library( precond SHARED ../../precond/asm_precond.cpp )
target_link_libraries( precond ngcomp ngsolve )
test_pde( asm.pde ${TEST_DIR}/precond ${TEST_DIR}/precond/libprecond.so )
