#include "Window.h"

void Window::init()
{
	if (glfwInit() != GLFW_TRUE)
	{
		ERROR("Can't init glfw.");
	}
	glfwSetErrorCallback(&Window::glfwError);
	glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
	glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);

	windowInit();
}

bool Window::windowShouldClose() const
{
	return glfwWindowShouldClose(_glfwWindow.get());
}

void Window::pollEvents()
{
	glfwPollEvents();
}

std::pair<const char**, std::uint32_t> Window::windowGetRequiredInstanceExtensions()
{
	std::uint32_t glfwExtensionCount = 0;
	const char** glfwExtensions;

	glfwExtensions = glfwGetRequiredInstanceExtensions(&glfwExtensionCount);
	return std::pair<const char**, std::uint32_t>(glfwExtensions, glfwExtensionCount);
}

void Window::windowInit()
{
	_glfwWindow.reset(glfwCreateWindow(800, 600, "cowboy-engine", nullptr, nullptr));
}

void Window::glfwError(int error, const char* description)
{
	WARNING(description);
}

Window::~Window()
{
	glfwTerminate();
}

void glfwDeleter::operator()(GLFWwindow* window)
{
	WARNING("GLFW Window deleted");
	glfwDestroyWindow(window);
}
