#include "Window.h"
#include "GLFW/glfw3.h"

void Window::init(WindowInfos windowInfos)
{
    _windowInfos = windowInfos;
	if (glfwInit() != GLFW_TRUE)
	{
		ERROR("Failed to initalize glfw.");
	}
	glfwSetErrorCallback(&Window::glfwError);

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
    //glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_FALSE);

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

const WindowInfos Window::windowInfos() const
{
    return _windowInfos;
}

void Window::swapBuffers()
{
    glfwSwapBuffers(_glfwWindow.get());
}

void Window::windowInit()
{
    _glfwWindow.reset(glfwCreateWindow(_windowInfos.x, _windowInfos.y, _windowInfos.title.c_str(), nullptr, nullptr));
	if (_glfwWindow == nullptr)
	{
        ERROR("Failed to create glfw window");		
	}
    glfwMakeContextCurrent(_glfwWindow.get());
}

void Window::glfwError([[maybe_unused]] int error, const char* description)
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
