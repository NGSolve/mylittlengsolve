project(mylittlengsolve)

cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)

set(CMAKE_INSTALL_PREFIX ${NGSOLVE_DIR} CACHE PATH "Install dir" FORCE)
set(PYTHON_PACKAGE_INSTALL_PATH ${NGSOLVE_INSTALL_DIR_PYTHON} CACHE PATH "Python install dir" FORCE)

add_library(myngsolve SHARED
  1_myFEM/myElement.cpp 1_myFEM/myFESpace.cpp 1_myFEM/myIntegrator.cpp
  1_myFEM/myPreconditioner.cpp
  2_myHOFEM/myHOFESpace.cpp 2_myHOFEM/myHOElement.cpp
  4_utility_functions/demo_coupling.cpp 4_utility_functions/myAssembling.cpp
  )

target_link_libraries(myngsolve ngsolve)
if(WIN32)
	target_compile_definitions(myngsolve PRIVATE MYNGS_DLL_HEADER=__declspec\(dllexport\))
else()
	target_compile_definitions(myngsolve PRIVATE MYNGS_DLL_HEADER)
endif()


add_ngsolve_python_module(myngspy myngspy.cpp)
set_target_properties(myngspy PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/${NETGEN_PYTHON_RPATH}")
target_link_libraries(myngspy myngsolve)

install(TARGETS myngspy DESTINATION ${PYTHON_PACKAGE_INSTALL_PATH})
