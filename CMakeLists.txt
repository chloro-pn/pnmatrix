cmake_minimum_required(VERSION 3.1)

project(pnmatrix C CXX)

enable_testing()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()

option (MAKE_TESTS "if you want to make tests" ON)
if(MAKE_TESTS)
	add_subdirectory("${PROJECT_SOURCE_DIR}/test")
endif()

option (MAKE_EXAMPLES "if you want to make examples" ON)
if(MAKE_EXAMPLES)
	add_subdirectory("${PROJECT_SOURCE_DIR}/examples")
endif()
