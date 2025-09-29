cmake_minimum_required(VERSION 3.5)
project(numeric_methods)

set(CMAKE_CXX_STANDARD 23)

include_directories(lab1/utils)

file(GLOB MATRIX_SRC "${CMAKE_SOURCE_DIR}/utils/*.c")

add_subdirectory(lab1)
