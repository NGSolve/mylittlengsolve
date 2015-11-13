add_library( mixed SHARED mixed_diffusion.cpp hybridmixed.cpp )
target_link_libraries( mixed ngcomp ngsolve )

install( TARGETS mixed DESTINATION lib )
install( FILES hybridmixed.pde mixed.pde square.in2d square.vol
         DESTINATION share/mylittlengsolve/mixed
       )
