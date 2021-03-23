#pragma once

#include "Shader.h"
#include "Window.h"
#include "Components.h"

class Renderer
{
public:
	void init(std::shared_ptr<Window> window);
    void prePass();
    void initMesh(Mesh& mesh);
    void drawMesh(Mesh& mesh);
    void freeMesh(Mesh& mesh);

private:
	std::shared_ptr<Window> _window = nullptr;

    Shader* shaderProgram;
};

