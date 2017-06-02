add_library(mixed SHARED mixed_diffusion hybridmixed.cpp)
target_link_libraries(mixed solve)
