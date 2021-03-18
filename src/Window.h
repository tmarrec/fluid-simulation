#pragma once
#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>
#include <memory>
#include <iostream>

#include "utils.h"

struct glfwDeleter
{
	void operator()(GLFWwindow* window);
};

class Window
{
public:
	void init();
	bool windowShouldClose() const;
	void pollEvents();
	std::pair<const char**, std::uint32_t> windowGetRequiredInstanceExtensions();
	~Window();



private:
	void windowInit();
	static void glfwError(int error, const char* description);

	std::unique_ptr<GLFWwindow, glfwDeleter> _glfwWindow = nullptr;
};
