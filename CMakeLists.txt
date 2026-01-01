cmake_minimum_required(VERSION 3.20)
project(numeric_methods)

set(CMAKE_CXX_STANDARD 23)

include_directories(utils)

file(GLOB MATRIX_SRC "${CMAKE_SOURCE_DIR}/utils/*.cpp")

add_subdirectory(lab1)
add_subdirectory(lab2)
add_subdirectory(lab3)
add_subdirectory(lab4)
