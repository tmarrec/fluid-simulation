#include "Physics.h"

void Physics::Update(float dt)
{
    for (auto const& entity : mEntities)
    {
        auto& transform = gCoordinator.GetComponent<Transform>(entity);
        transform.position += glm::vec3(0.1f, 0, 0);
    }
}
