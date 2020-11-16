#include <iostream>

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>

#include "ecs/Coordinator.h"
#include "utils.h"

Coordinator ECS;

int main()
{
	PRINT_TITLE();

	ECS.Init();
	std::cout << ":)" << std::endl;

	glfwInit();

    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    GLFWwindow* window = glfwCreateWindow(800, 600, "Vulkan window", nullptr, nullptr);

    uint32_t extensionCount = 0;
    vkEnumerateInstanceExtensionProperties(nullptr, &extensionCount, nullptr);

    std::cout << extensionCount << " extensions supported\n";
	glm::mat4 matrix;
    glm::vec4 vec;
    auto test = matrix * vec;


    while(!glfwWindowShouldClose(window)) {
        glfwPollEvents();
    }

    glfwDestroyWindow(window);

    glfwTerminate();

	return EXIT_SUCCESS;
}