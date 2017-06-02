add_library(precond SHARED asm_precond.cpp)
target_link_libraries(precond solve)
