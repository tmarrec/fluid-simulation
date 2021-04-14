#pragma once
#include <array>
#include <cstdint>
#include <vector>

#include "glm/fwd.hpp"
#include "glm/gtc/constants.hpp"
#include "utils.h"
#include "types.h"
#include "Shader.h"

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
    std::uint32_t TBO = 0;

    RenderMode renderMode = TRIANGLES;
};

struct Camera
{
    float yaw;
    float pitch;
    float speed;
    float FOV;
    glm::vec3 front = {0.0f, 0.0f, 1.0f};
    glm::vec3 up = {0.0f, 1.0f, 0.0f};
    glm::mat4 projection = glm::zero<glm::mat4>();
    Transform transform;
};

struct Material
{
    Shader shader;
    bool hasTexture = false;
    std::uint32_t texture = 0;
    std::uint32_t TBO = 0;
    std::vector<float> texCoords = {};
};

struct Fluid2D
{
    float viscosity;
    float dt;
    Entity entity;
    std::uint32_t N = 64;

    std::vector<float> velocityFieldX = {};
    std::vector<float> velocityFieldY = {};

    std::vector<float> velocityFieldPrevX = {};
    std::vector<float> velocityFieldPrevY = {};

    std::vector<float> substanceField = {};
    std::vector<float> substanceFieldPrev = {};

    std::uint32_t IX(const std::uint32_t x, int y) const
    { 
        return x + y * (N+2);
    };
};

struct Fluid3D
{
    float viscosity;
    float dt;
    Entity entity;
    std::uint32_t N = 16;

    std::vector<float> velocityFieldX = {};
    std::vector<float> velocityFieldY = {};
    std::vector<float> velocityFieldZ = {};

    std::vector<float> velocityFieldPrevX = {};
    std::vector<float> velocityFieldPrevY = {};
    std::vector<float> velocityFieldPrevZ = {};

    std::vector<float> substanceField = {};
    std::vector<float> substanceFieldPrev = {};

    std::uint32_t IX(const std::uint32_t x, const std::uint32_t y, const std::uint32_t z) const
    { 
        return x + y * (N+2) + z * ((N+2)*(N+2));
    };
};
