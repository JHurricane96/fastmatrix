cmake_minimum_required(VERSION 3.10.1)
project(fastmatrix)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

set(CMAKE_CXX_STANDARD 17)

option(BUILD_TESTS "Set to ON to build tests" ON)
option(BUILD_BENCHMARKS "Set to ON to build benchmarks" ON)

set(FASTMATRIX_PATH .)

if(BUILD_TESTS)
  add_executable(tests tests.cpp)
endif()

if(BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif()
