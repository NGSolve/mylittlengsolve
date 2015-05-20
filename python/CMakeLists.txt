add_library( myngspy SHARED myngspy.cpp )
target_link_libraries( myngspy ngfem ngcomp ngsolve )

install( TARGETS myngspy DESTINATION lib )
install( FILES myngspy.pde square.vol test.py
         DESTINATION share/mylittlengsolve/python
       )
