

cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)

add_library(mixed SHARED mixed_diffusion.cpp hybridmixed.cpp)
target_link_libraries(mixed ngsolve)
