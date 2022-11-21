cmake_minimum_required(VERSION 3.23)
project(numerical_methods)

set(CMAKE_CXX_STANDARD 14)

add_executable(numerical_methods main.cpp mh.cpp mh.h)
