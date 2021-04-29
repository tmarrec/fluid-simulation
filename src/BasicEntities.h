#pragma once

#include "ecs/Coordinator.h"
#include "Components.h"
#include "Shader.h"
#include "Renderer.h"

extern Coordinator gCoordinator;

class BasicEntities
{
public:
    static void initBasicEntities(std::shared_ptr<Renderer> renderer);
    static void addPlane(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
    static void addCube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
    static void addLineCube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
    static Entity addVector(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
    static void addFluid3D(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);

private:
    static void addTransform(const Entity& entity, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
    static Mesh _vector;
    static Mesh _plane;
    static Mesh _cube;
    static std::shared_ptr<Renderer> _renderer;
};
