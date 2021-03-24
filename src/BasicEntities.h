#pragma once

#include "ecs/Coordinator.h"
#include "Components.h"
#include "Shader.h"

extern Coordinator gCoordinator;

class BasicEntities
{
public:
    static void addCube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
    static void addVector(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);

private:
    static void addTransform(Entity& entity, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
};
