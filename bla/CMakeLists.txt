cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS /usr $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)

add_executable(matvec matvec.cpp) 
target_link_libraries(matvec ngbla) 

add_executable(timing timing.cpp) 
target_link_libraries(timing ngbla) 
