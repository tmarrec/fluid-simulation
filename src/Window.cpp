#include "Window.h"

void Window::init()
{
    if (glfwInit() != GLFW_TRUE)
    {
        ERROR("Failed to initalize glfw.");
    }
    glfwSetErrorCallback(&Window::glfwError);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_FALSE);

    windowInit();

    glfwSetKeyCallback(_glfwWindow.get(), Input::keyCallback);
    glfwSetCursorPosCallback(_glfwWindow.get(), Input::cursorPositionCallback);
    if (glfwRawMouseMotionSupported())
    {
        glfwSetInputMode(_glfwWindow.get(), GLFW_RAW_MOUSE_MOTION, GLFW_TRUE);
    }
}

bool Window::windowShouldClose() const
{
    return glfwWindowShouldClose(_glfwWindow.get());
}

void Window::pollEvents()
{
    glfwPollEvents();
}

void Window::swapBuffers()
{
    glfwSwapBuffers(_glfwWindow.get());
}

void Window::windowInit()
{
    _glfwWindow.reset(glfwCreateWindow(Config::width, Config::height,
                "fluid-simulation - Tristan Marrec 2021", nullptr, nullptr));
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

