add_executable( demo_std  demo_std.cpp )
target_link_libraries( demo_std ngstd )

add_executable( demo_bla demo_bla.cpp )
target_link_libraries( demo_bla ngstd ngbla )

add_executable( demo_fem demo_fem.cpp )
target_link_libraries( demo_fem ngfem ngstd )

add_executable( demo_comp demo_comp.cpp )
target_link_libraries( demo_comp ngcomp ngsolve ngla ngfem  ngstd nglib interface)

add_executable( demo_comp1d demo_comp1d.cpp netgen_dummy.cpp )
target_link_libraries( demo_comp1d ngcomp ngsolve ngla ngfem  ngstd nglib interface)

add_executable( demo_solve demo_solve.cpp )
target_link_libraries( demo_solve ngcomp ngsolve ngla ngfem  ngstd nglib interface)

add_executable( demo_parallel demo_parallel.cpp )
target_link_libraries( demo_parallel ngstd )

install( TARGETS 
         demo_std demo_bla demo_fem demo_comp demo_comp1d demo_solve demo_parallel
         DESTINATION share/mylittlengsolve/standalone
       )
