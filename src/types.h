#pragma once

#include <bits/c++config.h>
#include <bitset>
#include <cstdint>
#include <limits>

// ECS
using Entity = std::uint16_t;
const Entity MAX_ENTITIES = std::numeric_limits<Entity>::max();
using ComponentType = std::uint8_t;
const ComponentType MAX_COMPONENTS = std::numeric_limits<ComponentType>::max();
using Signature = std::bitset<MAX_COMPONENTS>;

struct WindowInfos
{
    std::string title;
    int x;
    int y;
};
