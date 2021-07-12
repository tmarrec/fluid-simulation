#pragma once
#include <array>
#include <cstdint>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include "glm/fwd.hpp"
#include "glm/gtc/constants.hpp"
#include "utils.h"
#include "types.h"
#include "Shader.h"

enum Solver
{
    CG,
    PCG
};

enum Advection
{
    SEMI_LAGRANGIAN,
    MACCORMACK
};

struct Transform
{
    glm::vec3 position;
    glm::vec3 rotation;
    glm::vec3 scale;
};

struct Mesh
{
    bool initialized = false;
    bool is2D = false;

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
    Shader shader = {};
    bool hasTexture = false;
    bool is2D = false;
    bool noShader = false;
    float absorption = 100.0f;
    glm::vec3 lightIntensity = glm::vec3(1.0, 1.0, 1.0);
    std::uint32_t texture = 0;
    std::uint32_t TBO = 0;
    std::vector<float> texCoords = {};
};


struct Fluid3D
{
    Entity entity;
};
