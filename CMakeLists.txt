cmake_minimum_required(VERSION 3.6)
project(untitled)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -lgomp") # uncomment to activate OPENMP

set(SOURCE_FILES main.cpp somClasses.cpp globals.h somTimer.cpp)
add_executable(untitled ${SOURCE_FILES})