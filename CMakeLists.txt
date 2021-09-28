CMAKE_MINIMUM_REQUIRED (VERSION 3.18.4)

PROJECT (fluid-simulation VERSION "0.1")

SET (ALLOW_DUPLICATE_CUSTOM_TARGETS TRUE)

OPTION (DISABLE_GUI "Disable window" OFF)
IF (DISABLE_GUI)
    ADD_DEFINITIONS(-DDISABLE_GUI)
    MESSAGE ("[fluid-simulation] GUI is disabled")
ENDIF()

# ╔═════════════════════════╗
# ║ C++ Compilation options ║
# ╚═════════════════════════╝

SET (CMAKE_CXX_STANDARD_REQUIRED ON)
SET (CMAKE_CXX_STANDARD 20)
SET (CMAKE_CXX_FLAGS "-std=c++20 -Wall -pedantic")
SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g -DDEBUG -lgomp -fopenmp")
SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -lgomp -fopenmp")


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

    # Main
 	src/main.cpp

    # Fluids
    src/Fluids.h
    src/Fluids.cpp
    src/ConjugateGradient.h
    src/ConjugateGradient.cpp
    src/StaggeredGrid.h

    src/Advect.h
    src/Advect.cpp
    src/Project.h
    src/Project.cpp
    src/MarchingCube.h
    src/MarchingCube.cpp

    # Engine
    src/Simulation.h
	src/Renderer.h
    src/Shader.h
	src/Window.h
    src/Input.h

    src/Simulation.cpp
	src/Renderer.cpp
    src/Shader.cpp
	src/Window.cpp
    src/Input.cpp

    # Export images
    extern/stb/stb_image_write.h

    # Read config file
    src/config.h
    src/config.cpp
    extern/inipp/inipp/inipp.h

    # Export .ply
    extern/tinyply/source/tinyply.h
)

# ╔═════════╗
# ║ Linking ║
# ╚═════════╝

TARGET_LINK_LIBRARIES (${PROJECT_NAME} OpenGL::GL)
TARGET_LINK_LIBRARIES (${PROJECT_NAME} glfw)
TARGET_LINK_LIBRARIES (${PROJECT_NAME} glad)

TARGET_INCLUDE_DIRECTORIES (${PROJECT_NAME} SYSTEM PRIVATE "extern/glm")
TARGET_INCLUDE_DIRECTORIES (${PROJECT_NAME} SYSTEM PRIVATE "extern/stb")
TARGET_INCLUDE_DIRECTORIES (${PROJECT_NAME} SYSTEM PRIVATE "extern/inipp/inipp")
TARGET_INCLUDE_DIRECTORIES (${PROJECT_NAME} SYSTEM PRIVATE "extern/tinyply/source")

FILE (COPY src/shaders DESTINATION .)
FILE (COPY src/config.ini DESTINATION .)

