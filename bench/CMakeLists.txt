cmake_minimum_required(VERSION 3.9.6)
project(fastmatrix)

add_executable(bench bench.cpp)
target_link_libraries(bench rt)
target_include_directories(bench PRIVATE
  ${CMAKE_SOURCE_DIR}/${FASTMATRIX_PATH})
