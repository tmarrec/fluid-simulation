#pragma once
#include <GLFW/glfw3.h>
#include <memory>
#include <iostream>

#include "types.h"
#include "utils.h"
#include "Input.h"

#ifdef DEBUG_GUI
#include "imgui_impl_glfw.h"
#endif

struct glfwDeleter
{
	void operator()(GLFWwindow* window);
};

class Window
{
public:
	void init(WindowInfos windowInfos);
	bool windowShouldClose() const;
	void pollEvents();
    void swapBuffers();
    const WindowInfos windowInfos() const;
	~Window();

#ifdef DEBUG_GUI
    void setupImgui() const;
#endif

private:
	void windowInit();
	static void glfwError(int error, const char* description);

	std::unique_ptr<GLFWwindow, glfwDeleter> _glfwWindow = nullptr;
    WindowInfos _windowInfos;
};
