project(mylittlengsolve)

cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)


add_ngsolve_python_module(myngspy myngspy.cpp
  1_myFEM/myElement.cpp 1_myFEM/myFESpace.cpp
  1_myFEM/myPreconditioner.cpp 1_myFEM/myCoefficient.cpp
  1_myFEM/myIntegrator.cpp
  2_myHOFEM/myHOFESpace.cpp 2_myHOFEM/myHOElement.cpp
  4_utility_functions/demo_coupling.cpp 4_utility_functions/myAssembling.cpp
  )

# check if CMAKE_INSTALL_PREFIX is set by user, otherwise install to user-local python modules directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import site; print(site.USER_SITE)" OUTPUT_VARIABLE install_dir OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_INSTALL_PREFIX ${install_dir} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

message("The python package will be installed to: ${CMAKE_INSTALL_PREFIX}")
install(TARGETS myngspy DESTINATION .)
