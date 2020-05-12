cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0015 NEW)

set (CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -fopenmp -Wextra")

execute_process(COMMAND ln -s ../shaders)

set(appName hecate)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLM REQUIRED)
find_package(Boost COMPONENTS filesystem REQUIRED )

include_directories(${OPENGL_INCLUDE_DIRS})
include_directories(${GLUT_INCLUDE_DIRS})
include_directories(${GLM_INCLUDE_DIRS})
include_directories("${CMAKE_SOURCE_DIR}/gl3w")
include_directories( ${Boost_INCLUDE_DIR} )

link_directories(${OPENGL_LIBRARY_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})

link_libraries(dl)

set(imgui_files imgui/imconfig.h imgui/imgui.cpp imgui/imgui_demo.cpp imgui/imgui_draw.cpp imgui/imgui.h imgui/imgui_internal.h imgui/imgui_widgets.cpp imgui/imstb_rectpack.h imgui/imstb_textedit.h imgui/imstb_truetype.h imgui/imgui_impl_opengl3.h imgui/imgui_impl_opengl3.cpp imgui/imgui_impl_glut.h imgui/imgui_impl_glut.cpp)

set(gl3w_files gl3w/gl3w.c gl3w/GL/gl3w.h)


file(GLOB hecate_SRC
    "*.h"
    "*.cpp"
    "geometry/*.cpp"
    "geometry/*.h"
    "lodepng/lodepng.h"
    "lodepng/lodepng.cpp"
)

add_executable(${appName} ${imgui_files} ${gl3w_files} ${hecate_SRC})

target_link_libraries(${appName} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${Boost_LIBRARIES})



