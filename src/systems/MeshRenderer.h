#pragma once
#include <cstdint>
#include <memory>

#include "../ecs/Coordinator.h"
#include "../Components.h"
#include "../Renderer.h"

extern Coordinator gCoordinator;

class MeshRenderer : public System
{
public:
    void init(std::shared_ptr<Renderer> renderer);
    void update();
    ~MeshRenderer();

private:
    std::shared_ptr<Renderer> _renderer = nullptr;
};
