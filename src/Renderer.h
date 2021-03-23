#pragma once

#include "Shader.h"
#include "Window.h"
#include "Components.h"
#include "utils.h"

class Renderer
{
public:
	void init(std::shared_ptr<Window> window);
    void prePass();
    void initMesh(Mesh& mesh);
    void drawMesh(Mesh& mesh);
    void freeMesh(Mesh& mesh);
    void useShader(Shader& shader, Camera& camera, Transform& transform);

private:
	std::shared_ptr<Window> _window = nullptr;
};

