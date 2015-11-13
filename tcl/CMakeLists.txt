add_library( tcldemo SHARED tcldemo.cpp )
target_link_libraries( tcldemo ngsolve )

install( TARGETS tcldemo DESTINATION lib )
install( FILES square.vol tcldemo.cpp tcldemo.pde tcldemo.tcl
         DESTINATION share/mylittlengsolve/tcl
       )
