#include "MeshRenderer.h"
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

MeshRenderer::~MeshRenderer()
{
    ASSERT(_renderer, "MeshRenderer needs to be initialized");
    for (auto const& entity : mEntities)
    {
        auto& mesh = gCoordinator.GetComponent<Mesh>(entity);
        if (mesh.initialized)
        {
            _renderer->freeMesh(mesh);
        }
    }
}

