#pragma once

#include <iostream>
#include <map>

#include "./GLFW/glfw3.h"

class Input
{
 public:
    static bool keyIsDown(int key);
    static void keyCallback(
            GLFWwindow* window,
            int key,
            int scancode,
            int action,
            int mob
        );
    static void cursorPositionCallback(
            GLFWwindow* window,
            double xpos,
            double ypos
        );
    static void updateMouseMovements();
    static float mouseOffsetX;
    static float mouseOffsetY;

 private:
    static std::map<int, bool> _keysStatus;
    static double _lastMouseX;
    static double _lastMouseY;
    static bool _focused;
};
