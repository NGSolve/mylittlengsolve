project(mylittlengsolve)

cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)

set(CMAKE_PREFIX_PATH ${NGSOLVE_INSTALL_DIR_PYTHON} CACHE PATH "Install dir")

add_ngsolve_python_module(myngspy myngspy.cpp
  1_myFEM/myElement.cpp 1_myFEM/myFESpace.cpp 1_myFEM/myIntegrator.cpp
  1_myFEM/myPreconditioner.cpp
  2_myHOFEM/myHOFESpace.cpp 2_myHOFEM/myHOElement.cpp
  4_utility_functions/demo_coupling.cpp 4_utility_functions/myAssembling.cpp
  )

target_link_libraries (myngspy ngsolve)

message("Package will be installed to: ${CMAKE_PREFIX_PATH}")
install(TARGETS myngspy DESTINATION ${CMAKE_PREFIX_PATH})
