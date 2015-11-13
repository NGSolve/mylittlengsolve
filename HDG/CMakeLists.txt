add_library( HDG SHARED HDGLaplace.cpp HDGElasticity.cpp )
target_link_libraries( HDG ngcomp ngsolve )

install( TARGETS HDG DESTINATION lib )
install( FILES
         HDGElasticity.pde HDGLaplace.pde pipe.geo pipe.vol.gz
         plate.geo plate.vol.gz square.in2d square.vol.gz
         DESTINATION share/mylittlengsolve/HDG
       )
