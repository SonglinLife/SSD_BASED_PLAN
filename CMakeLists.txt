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

add_executable(SSD_Based_Plan main.cpp)

set(TBB_TEST OFF CACHE BOOL "TBB Test" FORCE)
add_subdirectory(oneTBB)
target_link_libraries(SSD_Based_Plan PUBLIC tbb)