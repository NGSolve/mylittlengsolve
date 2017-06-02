add_library(HDG SHARED HDGLaplace.cpp HDGElasticity)
target_link_libraries(HDG solve)
