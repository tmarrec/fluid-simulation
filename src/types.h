#pragma once

#include <bitset>
#include <cstdint>

// ECS
using Entity = std::uint16_t;
const Entity MAX_ENTITIES = 65536;
using ComponentType = std::uint8_t;
const ComponentType MAX_COMPONENTS = 256;
using Signature = std::bitset<MAX_COMPONENTS>;