cmake_minimum_required(VERSION 2.8)
find_program( NGSCXX ngscxx )
if(NOT NGSCXX)
    message( FATAL_ERROR "Could not find ngscxx." )
endif(NOT NGSCXX)
set(CMAKE_CXX_COMPILER ${NGSCXX})
set(CMAKE_CXX_FLAGS "")

project(MyLittleNGSolve)

enable_testing()
include(CTest)

set( TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/original_files )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TEST_DIR} )
file(COPY ${PROJECT_SOURCE_DIR}/../ DESTINATION ${TEST_DIR} )

add_library( myngsolve SHARED ../all_in_one.cpp ../demo_instat.cpp ../demo_stokes.cpp ../myElement.cpp	 
            ../myHOElement.cpp ../myIntegrator.cpp ../demo_coupling.cpp ../demo_coupling_adv.cpp 
            ../demo_nonlinear.cpp ../myFESpace.cpp ../myHOFESpace.cpp ../myPreconditioner.cpp	 
            ../myAssembling.cpp ../linhypDG.cpp )

target_link_libraries( myngsolve ngfem ngcomp ngsolve )

function(test_pde pde_file working_dir )
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${pde_file}.appx")
        file(READ ${CMAKE_CURRENT_SOURCE_DIR}/${pde_file}.appx CONTENTS)
        file(APPEND ${working_dir}/${pde_file} ${CONTENTS})
    endif()      
    add_test(NAME ${pde_file}  COMMAND pde_tester "${working_dir}/${pde_file}" WORKING_DIRECTORY ${working_dir})
    set_tests_properties ( ${pde_file} PROPERTIES TIMEOUT 30 )
    set_tests_properties ( ${pde_file} PROPERTIES ATTACHED_FILES "${working_dir}/ng.prof" )
    set_tests_properties ( ${pde_file} PROPERTIES LABELS "standard")
    if(${ARGC} GREATER 2)
        set_tests_properties ( ${pde_file} PROPERTIES REQUIRED_FILES ${ARGN} )
    endif()
endfunction(test_pde)

set(pde_tests "all_in_one.pde;coupling.pde;linhyp.pde;my_assembling.pde;nonlinear.pde;simple.pde;stokes.pde" )
foreach(pde_file ${pde_tests})
  test_pde( ${pde_file} ${TEST_DIR} ${TEST_DIR}/libmyngsolve.so )
endforeach()
# instat.pde;precond.pde;  -> take too long

add_subdirectory(HDG)
add_subdirectory(mixed)
add_subdirectory(precond)

add_subdirectory(python)
add_subdirectory(standalone)
add_subdirectory(tcl)
