cmake_minimum_required(VERSION 3.10.1)
project(fastmatrix)

add_executable(bench bench.cpp)
target_link_libraries(bench rt)
target_include_directories(bench PRIVATE
  ${CMAKE_SOURCE_DIR}/${FASTMATRIX_PATH})
target_compile_options(bench PRIVATE "-O3;-DNDEBUG")
