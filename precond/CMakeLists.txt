add_library( precond SHARED asm_precond.cpp )
target_link_libraries( precond ngcomp ngsolve )

install( TARGETS precond DESTINATION lib )
install( FILES asm.pde square.in2d square.vol
         DESTINATION share/mylittlengsolve/precond
       )
