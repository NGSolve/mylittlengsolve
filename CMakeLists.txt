# project(mylittlengsolve)
# cmake_minimum_required(VERSION 3.1)

# find_package(NGSolve CONFIG REQUIRED
#   HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
# )

# set(CMAKE_INSTALL_PREFIX ${NGSOLVE_DIR} CACHE PATH "Install dir" FORCE)

# target_link_libraries(myngsolve ngsolve)
# if(WIN32)
# 	target_compile_definitions(myngsolve PRIVATE MYNGS_DLL_HEADER=__declspec\(dllexport\))
# else()
# 	target_compile_definitions(myngsolve PRIVATE MYNGS_DLL_HEADER)
# endif()

# add_subdirectory(1_Basic)
# add_subdirectory(2_Advanced)
# add_subdirectory(3_Legacy)

# install(TARGETS myngsolve RUNTIME DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} LIBRARY DESTINATION ${NGSOLVE_INSTALL_DIR_LIB} ARCHIVE DESTINATION ${NGSOLVE_INSTALL_DIR_LIB})
