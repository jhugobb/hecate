cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0015 NEW)

set (CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -fopenmp")

set(appName hec2png)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})

add_executable(${appName} main.cpp lodepng.h lodepng.cpp)

target_link_libraries(${appName})



