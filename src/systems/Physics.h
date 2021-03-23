#pragma once

#include "../ecs/Coordinator.h"
#include "../Components.h"

extern Coordinator gCoordinator;

class Physics : public System
{
public:
    void Update(float dt);
};
