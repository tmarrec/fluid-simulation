#pragma once

#include "ecs/Coordinator.h"
#include "Components.h"
#include "Shader.h"

extern Coordinator gCoordinator;

class BasicEntities
{
public:
    static void addPlane(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
    static void addCube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
    static void addLineCube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
    static void addVector(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
    static void addFluid2D(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);

private:
    static void addTransform(Entity& entity, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
};
