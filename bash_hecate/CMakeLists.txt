cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0015 NEW)

set (CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -fopenmp -Wextra")

set(appName hecate_bash)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../cmake_modules" ${CMAKE_MODULE_PATH})
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLM REQUIRED)

include_directories(${OPENGL_INCLUDE_DIRS})
include_directories(${GLUT_INCLUDE_DIRS})
include_directories(${GLM_INCLUDE_DIRS})
include_directories("${CMAKE_SOURCE_DIR}/../gl3w")

link_directories(${OPENGL_LIBRARY_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})

link_libraries(dl)

set(gl3w_files ../gl3w/gl3w.c ../gl3w/GL/gl3w.h)


file(GLOB hecate_SRC
    "../*.cpp"
    "../*.h"
    "*.h"
    "*.cpp"
    "../geometry/*.cpp"
    "../geometry/*.h"
    "../lodepng/lodepng.h"
    "../lodepng/lodepng.cpp"
)

list(REMOVE_ITEM hecate_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../main.cpp)
list(REMOVE_ITEM hecate_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../Scene.cpp)
list(REMOVE_ITEM hecate_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../Scene.h)
list(REMOVE_ITEM hecate_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../Application.cpp)
list(REMOVE_ITEM hecate_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../Application.h)

add_executable(${appName} ${gl3w_files} ${hecate_SRC})

target_link_libraries(${appName} ${GLUT_LIBRARIES})



