set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TEST_DIR}/python )
add_library( myngspy SHARED ../../python/myngspy.cpp )
target_link_libraries( myngspy ngfem ngcomp ngsolve )
test_pde( myngspy.pde ${TEST_DIR}/python ${TEST_DIR}/python/libmyngspy.so )

