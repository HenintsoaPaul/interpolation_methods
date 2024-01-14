cmake_minimum_required(VERSION 3.27)
project(interpolation C)

set(CMAKE_C_STANDARD 11)

add_executable(interpolation main.c)
target_link_libraries(interpolation m)