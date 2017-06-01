project(mylittlengsolve)
cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS /opt/netgen/lib/cmake /Applications/Netgen.app/Contents/Resource/CMake C:/netgen
)
        
add_library(myngsolve SHARED 
  all_in_one.cpp demo_instat.cpp demo_stokes.cpp myElement.cpp
  myHOElement.cpp myIntegrator.cpp demo_coupling.cpp demo_coupling_adv.cpp
  demo_nonlinear.cpp myFESpace.cpp myHOFESpace.cpp myPreconditioner.cpp
  myAssembling.cpp linhypDG.cpp periodic.cpp
  )

target_link_libraries(myngsolve solve)

install(TARGET myngsolve DESTINATION lib)
