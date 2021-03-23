#pragma once
#include <array>
#include <cstdint>
#include <vector>

#include "utils.h"

struct Transform
{
    glm::vec3 position;
    glm::vec3 rotation;
    glm::vec3 scale;
};

struct Mesh
{
    bool initialized = false;

    std::vector<float> vertices; 
    std::vector<float> normals; 
    std::vector<std::uint32_t> indices; 

    std::uint32_t VAO = 0;
    std::uint32_t VBO = 0;
    std::uint32_t NBO = 0;
    std::uint32_t EBO = 0;
};
