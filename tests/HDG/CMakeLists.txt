set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TEST_DIR}/HDG )
add_library( HDG SHARED ../../HDG/HDGLaplace.cpp ../../HDG/HDGElasticity.cpp )
target_link_libraries( HDG ngcomp ngsolve )
test_pde( HDGElasticity.pde ${TEST_DIR}/HDG libHDG.so )
test_pde( HDGLaplace.pde ${TEST_DIR}/HDG libHDG.so )

