cmake_minimum_required(VERSION 3.5)

project(cowboy-engine VERSION "0.0.4")

string(TIMESTAMP COMPILATION_DATE "%d/%m/%Y")


set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

cmake_policy(SET CMP0072 NEW)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -Wall -pedantic -Wextra -DGL_GLEXT_PROTOTYPES")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build)

find_package(Qt5 COMPONENTS OpenGL REQUIRED)
if (Qt5OpenGL_FOUND)
	message(STATUS "Found Qt ${Qt5OpenGL_VERSION}")
else (Qt5OpenGL_FOUND)
	message(EROR "Unable to find QTOpenGL library, you should probably change line 24 of CMakeLists.txt")
endif(Qt5OpenGL_FOUND)
find_package(OpenGL REQUIRED)

set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)

configure_file("src/config.h.in" "config.h")

add_executable(${PROJECT_NAME}
	src/config.h
 	src/main.cpp

	src/shapes/tiny_obj_loader.h

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

	src/Entity_Item.cpp
	src/Entity_Item.h

	src/Shape.cpp
	src/Shape.h

	src/Light.cpp
	src/Light.h

	src/Camera.cpp
	src/Camera.h

	src/shapes/BSpline.cpp
	src/shapes/BSpline.h

	src/shapes/Model.cpp
	src/shapes/Model.h
	src/shapes/Triangle.cpp
	src/shapes/Triangle.h
	src/shapes/Cube.cpp
	src/shapes/Cube.h
	src/shapes/UV_Sphere.cpp
	src/shapes/UV_Sphere.h
	src/shapes/Ico_Sphere.cpp
	src/shapes/Ico_Sphere.h
)

set(EXT_LIBS Qt5::OpenGL ${OPENGL_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${EXT_LIBS})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build)
file(COPY ${CMAKE_SOURCE_DIR}/src/shaders DESTINATION ${EXECUTABLE_OUTPUT_PATH})
file(COPY ${CMAKE_SOURCE_DIR}/src/models DESTINATION ${EXECUTABLE_OUTPUT_PATH})
file(COPY ${CMAKE_SOURCE_DIR}/src/images DESTINATION ${EXECUTABLE_OUTPUT_PATH})
