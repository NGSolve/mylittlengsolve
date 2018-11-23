

cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)

add_ngsolve_python_module(liblinhyp linhypDG.cpp)
# add_ngsolve_python_module(libequilibrate equilibrate_simple.cpp)
