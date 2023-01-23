cmake_minimum_required(VERSION 3.16)
add_executable(sedimentacao main.cc)

target_link_libraries(sedimentacao PRIVATE rcd)

find_package(OpenMP)
if (OPENMP_FOUND)
    target_link_libraries(sedimentacao PRIVATE OpenMP::OpenMP_CXX)
endif()

