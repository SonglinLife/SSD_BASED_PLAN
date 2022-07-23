cmake_minimum_required(VERSION 3.16)
project(SSD_Based_Plan)
set(CMAKE_CXX_STANDARD 17)

# OpenMP
find_package(OpenMP)

if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(FATAL_ERROR "no OpenMP supprot")
endif()

add_compile_options(-rdynamic)
add_executable(SSD_Based_Plan main.cpp)

target_include_directories(SSD_Based_Plan PUBLIC /data/wsl/SSD_Based_Plan/KaHIP/deploy)
target_link_directories(SSD_Based_Plan PUBLIC /data/wsl/SSD_Based_Plan/KaHIP/deploy)
target_link_directories(SSD_Based_Plan PUBLIC /usr/local/lib)
target_link_libraries(SSD_Based_Plan PUBLIC libkahip.a libmetis.so libGKlib.a)
target_compile_options(SSD_Based_Plan PUBLIC -fopenmp)

# add_compile_options(-fno-inline-small-functions)
