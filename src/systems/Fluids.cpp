#include "Fluids.h"

void Fluids::update()
{
    for (auto const& entity : mEntities)
    {
        auto& transform = gCoordinator.GetComponent<Transform>(entity);
    }
}
