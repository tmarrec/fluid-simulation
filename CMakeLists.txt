cmake_minimum_required(VERSION 3.11)

project(cowboy-engine VERSION "0.1.0")

string(TIMESTAMP COMPILATION_DATE "%d/%m/%Y")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#cmake_policy(SET CMP0072 NEW)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -Wall -pedantic -Wextra -pthread -DGL_GLEXT_PROTOTYPES")

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
set(CMAKE_AUTORCC on)

configure_file("src/config.h.in" "config.h")

add_executable(${PROJECT_NAME}
	src/config.h
 	src/main.cpp
	
	src/structs.h

	src/Shader.h
	src/Shader.cpp

	src/ECS.h
	src/DrawableComponent.h

	src/System.h
	src/System.cpp

	src/Message.h

	src/MessageBus.h
	src/MessageBus.cpp

	src/ui/MainWindow.h
	src/ui/MainWindow.cpp
	src/ui/MainWindow.ui

	src/Renderer.h
	src/Renderer.cpp

	src/ui/GlWidget.h
	src/ui/GlWidget.cpp
)

set(EXT_LIBS Qt5::OpenGL ${OPENGL_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${EXT_LIBS})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build)

file(COPY ${CMAKE_SOURCE_DIR}/src/shaders DESTINATION ${EXECUTABLE_OUTPUT_PATH})
#file(COPY ${CMAKE_SOURCE_DIR}/src/models DESTINATION ${EXECUTABLE_OUTPUT_PATH})
#file(COPY ${CMAKE_SOURCE_DIR}/src/images DESTINATION ${EXECUTABLE_OUTPUT_PATH})
