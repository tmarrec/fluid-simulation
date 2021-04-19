CMAKE_MINIMUM_REQUIRED (VERSION 3.18.4)

PROJECT (cowboy-engine VERSION "0.1")

# ╔═════════════════════════╗
# ║ C++ Compilation options ║
# ╚═════════════════════════╝

SET (CMAKE_CXX_STANDARD_REQUIRED ON)
SET (CMAKE_CXX_STANDARD 20)
SET (CMAKE_CXX_FLAGS "-std=c++20 -Wall -Wextra -pedantic")
SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -fopenmp")

# ╔════════════════════════════╗
# ║ OpenGL compilation options ║
# ╚════════════════════════════╝

FIND_PACKAGE (OpenGL REQUIRED)

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
	# Utils
	src/utils.h
	src/types.h

    # ECS
    src/ecs/System.h
    src/ecs/Coordinator.h
    src/ecs/EntityManager.h
    src/ecs/SystemManager.h
    src/ecs/ComponentArray.h
    src/ecs/ComponentManager.h

    # Components structs
    src/Components.h

    # Systems
    src/systems/Physics.h
    src/systems/Physics.cpp
    src/systems/MeshRenderer.h
    src/systems/MeshRenderer.cpp
    src/systems/Fluids.h
    src/systems/Fluids.cpp

    # Main
 	src/cowboy-engine.cpp

    # Engine
	src/Game.h
	src/Renderer.h
	src/Window.h
    src/Shader.h
    src/Input.h
    src/BasicEntities.h

	src/Game.cpp
	src/Renderer.cpp
	src/Window.cpp
    src/Shader.cpp
    src/Input.cpp
    src/BasicEntities.cpp
)

# ╔═════════╗
# ║ Linking ║
# ╚═════════╝

TARGET_LINK_LIBRARIES (${PROJECT_NAME} OpenGL::GL)
TARGET_LINK_LIBRARIES (${PROJECT_NAME} glfw)
TARGET_LINK_LIBRARIES (${PROJECT_NAME} glad)
TARGET_INCLUDE_DIRECTORIES (${PROJECT_NAME} PRIVATE "extern/glm")


FILE (COPY src/shaders DESTINATION .)

