CMAKE_MINIMUM_REQUIRED (VERSION 3.18.4)

PROJECT (cowboy-engine VERSION "0.1")

SET (ALLOW_DUPLICATE_CUSTOM_TARGETS TRUE)

OPTION (DEBUG_GUI "Debug interface" OFF)
IF (DEBUG_GUI)
    ADD_DEFINITIONS(-DDEBUG_GUI)
    MESSAGE ("[cowboy-engine] Using debug GUI")
ENDIF()

OPTION (DISABLE_GUI "Disable window" OFF)
IF (DISABLE_GUI)
    ADD_DEFINITIONS(-DDISABLE_GUI)
    MESSAGE ("[cowboy-engine] GUI is disabled")
ENDIF()

# ╔═════════════════════════╗
# ║ C++ Compilation options ║
# ╚═════════════════════════╝

SET (CMAKE_CXX_STANDARD_REQUIRED ON)
SET (CMAKE_CXX_STANDARD 20)
SET (CMAKE_CXX_FLAGS "-std=c++20 -Wall -Wextra -pedantic")
SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g -fopenmp")
SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -fopenmp")

# ╔════════════════════════════╗
# ║ OpenGL compilation options ║
# ╚════════════════════════════╝

FIND_PACKAGE (OpenGL REQUIRED)

# ╔════════════════════════════╗
# ║ eigen3 compilation options ║
# ╚════════════════════════════╝

INCLUDE_DIRECTORIES(SYSTEM extern/eigen)

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
    src/systems/MeshRenderer.h
    src/systems/MeshRenderer.cpp
    src/systems/Fluids.h
    src/systems/Fluids.cpp

    # Main
 	src/cowboy-engine.cpp

    # Engine
	src/Game.h
	src/Renderer.h
    src/Shader.h
    src/BasicEntities.h
	src/Window.h
    src/Input.h

	src/Game.cpp
	src/Renderer.cpp
    src/Shader.cpp
    src/BasicEntities.cpp
	src/Window.cpp
    src/Input.cpp

    # imgui
    extern/imgui/imgui.cpp
    extern/imgui/imgui_draw.cpp
    extern/imgui/imgui_tables.cpp
    extern/imgui/imgui_widgets.cpp
    extern/imgui/backends/imgui_impl_glfw.cpp
    extern/imgui/backends/imgui_impl_opengl3.cpp
    extern/implot/implot.cpp
    extern/implot/implot_items.cpp

    # stb
    extern/stb/stb_image_write.h
)

# ╔═════════╗
# ║ Linking ║
# ╚═════════╝

TARGET_LINK_LIBRARIES (${PROJECT_NAME} OpenGL::GL)
TARGET_LINK_LIBRARIES (${PROJECT_NAME} glfw)
TARGET_LINK_LIBRARIES (${PROJECT_NAME} glad)

TARGET_INCLUDE_DIRECTORIES (${PROJECT_NAME} SYSTEM PRIVATE "extern/imgui")
TARGET_INCLUDE_DIRECTORIES (${PROJECT_NAME} SYSTEM PRIVATE "extern/imgui/backends")
TARGET_INCLUDE_DIRECTORIES (${PROJECT_NAME} SYSTEM PRIVATE "extern/implot")
TARGET_INCLUDE_DIRECTORIES (${PROJECT_NAME} SYSTEM PRIVATE "extern/glm")
TARGET_INCLUDE_DIRECTORIES (${PROJECT_NAME} SYSTEM PRIVATE "extern/stb")

FILE (COPY src/shaders DESTINATION .)

