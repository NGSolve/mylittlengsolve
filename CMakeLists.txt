project(mylittlengsolve)
cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)

set(CMAKE_INSTALL_PREFIX ${NGSOLVE_DIR} CACHE PATH "Install dir" FORCE)
        
add_library(myngsolve SHARED 
  all_in_one.cpp demo_instat.cpp demo_stokes.cpp myElement.cpp
  myHOElement.cpp myIntegrator.cpp demo_coupling.cpp demo_coupling_adv.cpp
  demo_nonlinear.cpp myFESpace.cpp myHOFESpace.cpp myPreconditioner.cpp
  myAssembling.cpp linhypDG.cpp periodic.cpp
  )

target_link_libraries(myngsolve ngsolve)
if(WIN32)
	target_compile_definitions(myngsolve PRIVATE MYNGS_DLL_HEADER=__declspec\(dllexport\))
else()
	target_compile_definitions(myngsolve PRIVATE MYNGS_DLL_HEADER)
endif()

add_subdirectory(python)
# add_subdirectory(standalone)
add_subdirectory(mixed)
add_subdirectory(HDG)
add_subdirectory(precond)
add_subdirectory(equilibrate)
# add_subdirectory(tcl)

install(TARGETS myngsolve RUNTIME DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} LIBRARY DESTINATION ${NGSOLVE_INSTALL_DIR_LIB} ARCHIVE DESTINATION ${NGSOLVE_INSTALL_DIR_LIB})
