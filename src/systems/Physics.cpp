#include "Physics.h"

void Physics::update(float dt)
{
    for (auto const& entity : mEntities)
    {
        auto& transform = gCoordinator.GetComponent<Transform>(entity);
        transform.rotation += glm::vec3(20.0f, 15.0f, 10.0f) * dt;
    }
}
