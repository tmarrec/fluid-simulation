#include "MeshRenderer.h"
#include "GLFW/glfw3.h"
#include "glm/ext/matrix_transform.hpp"
#include "glm/trigonometric.hpp"

void MeshRenderer::init(std::shared_ptr<Renderer> renderer, Camera camera)
{
    _renderer = renderer;
    _camera = camera;
}

void MeshRenderer::update()
{
    ASSERT(_renderer, "MeshRenderer needs to be initialized");
    cameraMovements();
    for (auto const& entity : mEntities)
    {
        auto& transform = gCoordinator.GetComponent<Transform>(entity);
        auto& mesh = gCoordinator.GetComponent<Mesh>(entity);
        auto& material = gCoordinator.GetComponent<Material>(entity);

        if (!mesh.initialized)
        {
            _renderer->initMesh(mesh);
        }
        else
        {
            _renderer->useShader(material.shader, _camera, transform);
            _renderer->drawMesh(mesh);
        }
    }
}

void MeshRenderer::cameraMovements()
{
    if (KeyInput::isDown(GLFW_KEY_W))
    {
        _camera.transform.position.z += _camera.speed;
    }
    if (KeyInput::isDown(GLFW_KEY_A))
    {
        _camera.transform.position.x += _camera.speed;
    }
    if (KeyInput::isDown(GLFW_KEY_S))
    {
        _camera.transform.position.z -= _camera.speed;
    }
    if (KeyInput::isDown(GLFW_KEY_D))
    {
        _camera.transform.position.x -= _camera.speed;
    }
    if (KeyInput::isDown(GLFW_KEY_LEFT_SHIFT))
    {
        _camera.transform.position.y += _camera.speed;
    }
    if (KeyInput::isDown(GLFW_KEY_LEFT_CONTROL))
    {
        _camera.transform.position.y -= _camera.speed;
    }

}

MeshRenderer::~MeshRenderer()
{
    ASSERT(_renderer, "MeshRenderer needs to be initialized");
    for (const auto& entity : mEntities)
    {
        auto& mesh = gCoordinator.GetComponent<Mesh>(entity);
        if (mesh.initialized)
        {
            _renderer->freeMesh(mesh);
        }
    }
}

