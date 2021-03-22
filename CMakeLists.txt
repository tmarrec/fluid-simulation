CMAKE_MINIMUM_REQUIRED (VERSION 3.18.4)

PROJECT (cowboy-engine VERSION "0.1")

# ╔═════════════════════════╗
# ║ C++ Compilation options ║
# ╚═════════════════════════╝

SET (CMAKE_CXX_STANDARD_REQUIRED ON)
SET (CMAKE_CXX_STANDARD 17)
SET (CMAKE_CXX_FLAGS "-std=c++17 -Wall -Wextra -pedantic")
SET (CMAKE_CXX_FLAGS_DEBUG_INIT "${CMAKE_CXX_FLAGS} -O0 -g")
SET (CMAKE_CXX_FLAGS_RELEASE_INIT "${CMAKE_CXX_FLAGS} -O3")

# ╔══════════════════════════╗
# ║ glfw compilation options ║
# ╚══════════════════════════╝

SET (GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
SET (GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
SET (GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
ADD_SUBDIRECTORY (extern/glfw)

# ╔══════════════════════════╗
# ║ glad compilation options ║
# ╚══════════════════════════╝

ADD_SUBDIRECTORY (extern/glad)

# ╔═════════════╗
# ║ Vim options ║
# ╚═════════════╝

SET (CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ╔═════════╗
# ║ Sources ║
# ╚═════════╝

ADD_EXECUTABLE (${PROJECT_NAME}
	# Header only
	src/utils.h
	src/types.h

	# Source only
 	src/cowboy-engine.cpp

	# Header and source files
	src/Game.h
	src/Renderer.h
	src/Window.h

	src/Game.cpp
	src/Renderer.cpp
	src/Window.cpp
)

# ╔══════════════╗
# ║ glfw linking ║
# ╚══════════════╝

TARGET_LINK_LIBRARIES (${PROJECT_NAME} glfw)
FIND_PACKAGE (OpenGL REQUIRED)
TARGET_LINK_LIBRARIES (${PROJECT_NAME} OpenGL::GL)

# ╔══════════════╗
# ║ glad linking ║
# ╚══════════════╝

TARGET_LINK_LIBRARIES(${PROJECT_NAME} glad)
