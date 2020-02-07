cmake_minimum_required(VERSION 3.5)

project(cowboy-engine)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

cmake_policy(SET CMP0072 NEW)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic -Wextra -DGL_GLEXT_PROTOTYPES")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

find_package(Qt5 COMPONENTS OpenGL REQUIRED)
if (Qt5OpenGL_FOUND)
	message(STATUS "Found Qt ${Qt5OpenGL_VERSION}")
else (Qt5OpenGL_FOUND)
	message(EROR "Unable to find QTOpenGL library")
endif(Qt5OpenGL_FOUND)
find_package(OpenGL REQUIRED) # define OPENGL_LIBRARIES

add_executable(${PROJECT_NAME}
 	src/main.cpp
	
	src/MainWindow.cpp
	src/MainWindow.h

	src/GLWidget.h
	src/GLWidget.cpp

	src/OpenGL.cpp
	src/OpenGL.h

	src/Shader.cpp
	src/Shader.h
	
	src/ECS.cpp
	src/ECS.h

	src/Entity.cpp
	src/Entity.h

	src/shapes/Triangle.cpp
	src/shapes/Triangle.h
)

set(EXT_LIBS Qt5::OpenGL ${OPENGL_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${EXT_LIBS})
