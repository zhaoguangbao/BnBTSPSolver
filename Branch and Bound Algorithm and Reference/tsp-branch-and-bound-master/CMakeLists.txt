CMAKE_MINIMUM_REQUIRED(VERSION 3.4)
project(branch_and_bound)

include_directories(include)

set(SOURCES
    src/main.cpp
    src/Parser.cpp
    src/Instance.cpp
    src/Solver.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++14 -Wall -Wextra")

add_executable(branch_and_bound ${SOURCES})
