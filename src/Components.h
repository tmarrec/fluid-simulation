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
    Shader shader = {};
    bool hasTexture = false;
    bool noShader = false;
    float absorption = 64.0f;
    glm::vec3 lightIntensity = glm::vec3(1, 1, 1);
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
    Entity entity;
    float viscosity = 0.2f;
    float diffusion = 0;
    float dt = 0.0001f;
    std::uint32_t N = 12;
    std::uint8_t solveNb = 8;

    std::vector<float> velocityFieldX = {};
    std::vector<float> velocityFieldY = {};
    std::vector<float> velocityFieldZ = {};

    std::vector<float> velocityFieldPrevX = {};
    std::vector<float> velocityFieldPrevY = {};
    std::vector<float> velocityFieldPrevZ = {};

    std::vector<float> substanceField = {};
    std::vector<float> substanceFieldPrev = {};

    //Eigen::SparseMatrix<double> laplacian = Eigen::SparseMatrix<double>(N*N*N,N*N*N);

    std::uint32_t IX(const std::uint32_t x, const std::uint32_t y, const std::uint32_t z) const
    { 
        return x + y * (N+2) + z * ((N+2)*(N+2));
    };

    void init()
    {
        velocityFieldX.reserve((N+2)*(N+2)*(N+2));
        for (std::uint32_t i = 0; i < (N+2)*(N+2)*(N+2); ++i)
        {
            velocityFieldX.emplace_back(0);
        }
        velocityFieldPrevX = velocityFieldX;
        velocityFieldY = velocityFieldX;
        velocityFieldPrevY = velocityFieldX;
        velocityFieldZ = velocityFieldX;
        velocityFieldPrevZ = velocityFieldX;
        substanceField = velocityFieldX;
        substanceFieldPrev = velocityFieldX;
    };
};
