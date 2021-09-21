#include "Input.h"

std::map<int, bool> Input::_keysStatus;
double Input::_lastMouseX;
double Input::_lastMouseY;
float Input::mouseOffsetX;
float Input::mouseOffsetY;
bool Input::_focused = false;

// Key callback, mainly used to know if the windows is focused or not
void Input::keyCallback(
        GLFWwindow* window,
        int key,
        [[maybe_unused]] int scancode,
        int action,
        [[maybe_unused]] int mobs
    )
{
    switch(action)
    {
        case GLFW_PRESS:
            _keysStatus[key] = true;
            break;
        case GLFW_RELEASE:
            _keysStatus[key] = false;
            break;
    }
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        _focused = false;
    }
    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
    {
        if (!_focused)
        {
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
            _focused = true;
        }
        else
        {
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
            _focused = false;
        }
    }
}

// Key pressed, add it to the _keysStatus map
bool Input::keyIsDown(int key)
{
    bool result = false;
    std::map<int, bool>::iterator it = _keysStatus.find(key);
    if (it != _keysStatus.end())
    {
        result = _keysStatus[key];
    }
    return result;
}

// Get the mouse position
void Input::cursorPositionCallback(
        [[maybe_unused]] GLFWwindow* window,
        double xpos,
        double ypos
    )
{
    if (_focused)
    {
        mouseOffsetX = xpos-_lastMouseX;
        mouseOffsetY = ypos-_lastMouseY;
    }
    _lastMouseX = xpos;
    _lastMouseY = ypos;
}

// Reset mouse position offsets
void Input::updateMouseMovements()
{
    mouseOffsetX = 0.0f;
    mouseOffsetY = 0.0f;
}
