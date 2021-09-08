#pragma once

#include "types.h"
#include "utils.h"
#include "Input.h"

#include <GLFW/glfw3.h>
#include <memory>
#include <iostream>

struct glfwDeleter
{
	void operator()(GLFWwindow* window);
};

class Window
{
public:
	void init(const WindowInfos& windowInfos);
	bool windowShouldClose() const;
	void pollEvents();
    void swapBuffers();
    const WindowInfos windowInfos() const;
	~Window();

private:
	void windowInit();
	static void glfwError(int error, const char* description);

	std::unique_ptr<GLFWwindow, glfwDeleter> _glfwWindow = nullptr;
    WindowInfos _windowInfos {};
};
