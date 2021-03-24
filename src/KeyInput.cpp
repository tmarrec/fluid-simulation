#include "KeyInput.h"
#include "GLFW/glfw3.h"

std::map<int, bool> KeyInput::_keysStatus;

void KeyInput::keyCallback([[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action, [[maybe_unused]] int mobs)
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
    if (key == GLFW_KEY_W && action == GLFW_PRESS)
    {
        _keysStatus[key] = true;
    }
}

bool KeyInput::isDown(int key)
{
    bool result = false; 
    std::map<int, bool>::iterator it = _keysStatus.find(key);
    if (it != _keysStatus.end())
    {
        result = _keysStatus[key];
    }
    return result;
}
