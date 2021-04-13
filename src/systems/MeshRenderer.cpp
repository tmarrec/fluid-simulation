#include "MeshRenderer.h"

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
            _renderer->applyMaterial(material, _camera, transform);
            _renderer->drawMesh(mesh);
        }
    }
}

void MeshRenderer::cameraMovements()
{
    if (Input::keyIsDown(GLFW_KEY_W))
    {
        _camera.transform.position += _camera.front * _camera.speed;
    }
    if (Input::keyIsDown(GLFW_KEY_A))
    {
        _camera.transform.position -= glm::normalize(glm::cross(_camera.front, _camera.up)) * _camera.speed;
    }
    if (Input::keyIsDown(GLFW_KEY_S))
    {
        _camera.transform.position -= _camera.front * _camera.speed;
    }
    if (Input::keyIsDown(GLFW_KEY_D))
    {
        _camera.transform.position += glm::normalize(glm::cross(_camera.front, _camera.up)) * _camera.speed;
    }
    if (Input::keyIsDown(GLFW_KEY_LEFT_SHIFT))
    {
        _camera.transform.position.y += _camera.speed;
    }
    if (Input::keyIsDown(GLFW_KEY_LEFT_CONTROL))
    {
        _camera.transform.position.y -= _camera.speed;
    }
    cameraMouseMovements();
}

void MeshRenderer::cameraMouseMovements()
{
    float sensitivity = 1.0f;
    _camera.yaw += Input::mouseOffsetX*sensitivity;
    _camera.pitch = glm::clamp(_camera.pitch - Input::mouseOffsetY*sensitivity, -89.9f, 89.9f);

    glm::vec3 dir;
    dir.x = cos(glm::radians(_camera.yaw))*cos(glm::radians(_camera.pitch));
    dir.y = sin(glm::radians(_camera.pitch));
    dir.z = sin(glm::radians(_camera.yaw))*cos(glm::radians(_camera.pitch));
    dir = glm::normalize(dir);
    _camera.front = dir;

    Input::updateMouseMovements();
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

