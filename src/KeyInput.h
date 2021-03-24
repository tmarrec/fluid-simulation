#pragma once
#include "GLFW/glfw3.h"
#include <iostream>
#include <map>

class KeyInput
{
public:
    static bool isDown(int key);
    static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mobs);

private:
    static std::map<int, bool> _keysStatus;
};
