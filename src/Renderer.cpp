#include "Renderer.h"

void Renderer::init(std::shared_ptr<Window> window)
{
	_window = window;

    if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress)))
    {
        ERROR("Failed to initialize glad");
    }

    const WindowInfos windowInfos = _window->windowInfos();
    glViewport(0, 0, windowInfos.x, windowInfos.y);
}

void Renderer::pass()
{
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
}

Renderer::~Renderer()
{

}

