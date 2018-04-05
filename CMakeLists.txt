cmake_minimum_required (VERSION 2.6)

project(main)

# make sure that the default is a RELEASE

file(GLOB SOURCE "./*.h" "./*.cpp")
set (CMAKE_CXX_FLAGS "-O3 -std=c++11 -g")
#set (CMAKE_C_FLAGS "-std=c99")

add_executable (main ${SOURCE})

