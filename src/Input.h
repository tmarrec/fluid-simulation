#pragma once
#include "GLFW/glfw3.h"
#include <iostream>
#include <map>

class Input
{
public:
    static bool keyIsDown(int key);
    static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mobs);
    static void cursorPositionCallback(GLFWwindow* window, double xpos, double ypos);
    static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mobs); 
    static void updateMouseMovements();
    static float mouseOffsetX;
    static float mouseOffsetY;

private:
    static std::map<int, bool> _keysStatus;
    static double _lastMouseX;
    static double _lastMouseY;
    static bool _focused;
};
