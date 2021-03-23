#include "MeshRenderer.h"

void MeshRenderer::init(std::shared_ptr<Renderer> renderer)
{
    _renderer = renderer;
}

void MeshRenderer::update()
{
    ASSERT(_renderer, "MeshRenderer needs to be initialized");
    for (auto const& entity : mEntities)
    {
        auto& transform = gCoordinator.GetComponent<Transform>(entity);
        auto& mesh = gCoordinator.GetComponent<Mesh>(entity);

        if (!mesh.initialized)
        {
            _renderer->initMesh(mesh);
        }
        else
        {
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

