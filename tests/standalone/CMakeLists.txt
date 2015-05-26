set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR}/standalone )
add_executable( demo_std  ../../standalone/demo_std.cpp )
target_link_libraries( demo_std ngstd )

add_executable( demo_bla ../../standalone/demo_bla.cpp )
target_link_libraries( demo_bla ngstd ngbla )

add_executable( demo_fem ../../standalone/demo_fem.cpp )
target_link_libraries( demo_fem ngfem ngstd )

add_executable( demo_comp ../../standalone/demo_comp.cpp )
target_link_libraries( demo_comp ngcomp ngsolve ngla ngfem  ngstd nglib interface)

add_executable( demo_comp1d ../../standalone/demo_comp1d.cpp ../../standalone/netgen_dummy.cpp )
target_link_libraries( demo_comp1d ngcomp ngsolve ngla ngfem  ngstd nglib interface)

add_executable( demo_solve ../../standalone/demo_solve.cpp )
target_link_libraries( demo_solve ngcomp ngsolve ngla ngfem  ngstd nglib interface)

add_executable( demo_parallel ../../standalone/demo_parallel.cpp )
target_link_libraries( demo_parallel ngstd )
