cmake_minimum_required(VERSION 3.22)
project(untitled)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lgomp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -lgomp")

set(SOURCE_FILES main.cpp somClasses.cpp globals.h somTimer.cpp metrics.h)
add_executable(untitled ${SOURCE_FILES})