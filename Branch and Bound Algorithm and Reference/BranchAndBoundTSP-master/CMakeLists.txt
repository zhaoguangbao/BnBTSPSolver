cmake_minimum_required(VERSION 3.6)
project(BranchAndBoundTSP)

set(CMAKE_CXX_STANDARD 14)

file(GLOB_RECURSE BranchAndBoundTSP_SOURCES "src/*.cpp")
file(GLOB_RECURSE BranchAndBoundTSP_HEADERS "header/*.hpp")

set (BranchAndBoundTSP_INCLUDE_DIRS "")
foreach (_headerFile ${BranchAndBoundTSP_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND BranchAndBoundTSP_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES BranchAndBoundTSP_INCLUDE_DIRS)

set(GCC_COVERAGE_COMPILE_FLAGS "-std=c++14 -Wall -Wshadow  -Wextra -pedantic -g   -march=native  -Werror  -O2") #

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )
include_directories(${BranchAndBoundTSP_INCLUDE_DIRS})

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygenConfig.cfg ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
            ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif(DOXYGEN_FOUND)

add_executable(BranchAndBoundTSP ${BranchAndBoundTSP_SOURCES})

#set(SOURCE_FILES src/main.cpp)
#add_executable(BranchAndBoundTSP ${SOURCE_FILES})