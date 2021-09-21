#pragma once

#include <memory>
#include <iostream>

#include "./types.h"
#include "./utils.h"
#include "./config.h"
#include "./Input.h"
#include "./GLFW/glfw3.h"

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
    void swapBuffers();
    ~Window();

 private:
    void windowInit();
    static void glfwError(int error, const char* description);

    std::unique_ptr<GLFWwindow, glfwDeleter> _glfwWindow = nullptr;
};
