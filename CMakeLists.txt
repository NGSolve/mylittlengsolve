cmake_minimum_required(VERSION 2.8)

find_program( NGSCXX ngscxx )
if(NOT NGSCXX)
    message( FATAL_ERROR "Could not find ngscxx." )
endif(NOT NGSCXX)
set(CMAKE_CXX_COMPILER ${NGSCXX})
set(CMAKE_CXX_FLAGS "")
message("found ngscxx: ${CMAKE_CXX_COMPILER_INIT}")

project(MyLittleNGSolve)

#######################################################################
# find netgen/ngsolve
if(NOT INSTALL_DIR)
    get_filename_component(NETGEN_BINDIR ${NGSCXX} PATH)
    get_filename_component(NETGEN_BASEDIR ${NETGEN_BINDIR} PATH)
    set(INSTALL_DIR ${NETGEN_BASEDIR} CACHE PATH "Path to install MyLittleNGSolve (i.e. where netgen is installed)")
endif(NOT INSTALL_DIR)
set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

add_library( myngsolve SHARED all_in_one.cpp demo_instat.cpp demo_stokes.cpp myElement.cpp	 
            myHOElement.cpp myIntegrator.cpp demo_coupling.cpp demo_coupling_adv.cpp 
            demo_nonlinear.cpp myFESpace.cpp myHOFESpace.cpp myPreconditioner.cpp	 
            myAssembling.cpp linhypDG.cpp )


target_link_libraries( myngsolve ngfem ngcomp ngsolve )

install( TARGETS myngsolve DESTINATION lib )
install( FILES
         all_in_one.pde coupling.pde instat.pde linhyp.pde my_assembling.pde nonlinear.pde precond.pde
         simple.pde stokes.pde square.in2d square.vol
         DESTINATION share/mylittlengsolve
       )

add_subdirectory(HDG)
add_subdirectory(mixed)
add_subdirectory(precond)

# add_subdirectory(python)
# add_subdirectory(standalone)
# add_subdirectory(tcl)
