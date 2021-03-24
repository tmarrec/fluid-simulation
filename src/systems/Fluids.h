#pragma once

#include "../ecs/Coordinator.h"
#include "../Components.h"

extern Coordinator gCoordinator;

class Fluids : public System
{
public:
    void update();
};
