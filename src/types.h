#pragma once

#include <bitset>
#include <cstdint>

// ECS
using Entity = std::uint16_t;
const Entity MAX_ENTITIES = static_cast<Entity>(65536);
using ComponentType = std::uint8_t;
const ComponentType MAX_COMPONENTS = static_cast<ComponentType>(256);
using Signature = std::bitset<MAX_COMPONENTS>;

struct WindowInfos
{
    std::string title;
    int x;
    int y;
};
