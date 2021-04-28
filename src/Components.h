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
    float absorption = 60.0f;
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
    double viscosity = 3.33;
    double diffusion = 0.15;
    double dt = 0.00005;
    std::uint32_t N = 32;

    std::vector<double> velocityFieldX = {};
    std::vector<double> velocityFieldY = {};
    std::vector<double> velocityFieldZ = {};

    std::vector<double> velocityFieldPrevX = {};
    std::vector<double> velocityFieldPrevY = {};
    std::vector<double> velocityFieldPrevZ = {};

    std::vector<double> substanceField = {};
    std::vector<double> substanceFieldPrev = {};

    Eigen::SparseMatrix<double> laplacianProject;
    Eigen::SparseMatrix<double> laplacianViscosity;
    Eigen::SparseMatrix<double> laplacianDiffuse;

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

        double visc = dt * viscosity * N;
        double diff = dt * diffusion * N;

        laplacianProject = Eigen::SparseMatrix<double>(N*N*N, N*N*N);
        laplacianViscosity = Eigen::SparseMatrix<double>(N*N*N, N*N*N);
        laplacianDiffuse = Eigen::SparseMatrix<double>(N*N*N, N*N*N);

        laplacianProject.reserve(Eigen::VectorXi::Constant(N*N*N, 6));
        laplacianViscosity.reserve(Eigen::VectorXi::Constant(N*N*N, 6));
        laplacianDiffuse.reserve(Eigen::VectorXi::Constant(N*N*N, 6));

        std::vector<Eigen::Triplet<double>> tripletListProject;
        std::vector<Eigen::Triplet<double>> tripletListViscosity;
        std::vector<Eigen::Triplet<double>> tripletListDiffuse;

        tripletListProject.reserve(6*N*N*N);
        tripletListViscosity.reserve(6*N*N*N);
        tripletListDiffuse.reserve(6*N*N*N);

        for(std::uint32_t i = 0; i < N*N*N; ++i)
        {
            if (i > 0 && (i-1)%N != 0)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i-1, 1));
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i-1, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i-1, diff));
            }
            if (i >= N)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i-N, 1));
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i-N, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i-N, diff));
            }
            if (i >= N * N)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i-N*N, 1));
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i-N*N, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i-N*N, diff));
            }
            if (i < N*N*N-1 && (i+1)%N != 0)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i+1, 1));
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i+1, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i+1, diff));
            }
            if (i + N < N*N*N-1)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i+N, 1));
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i+N, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i+N, diff));
            }
            if (i + N*N < N*N*N-1)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i+N*N, 1));
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i+N*N, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i+N*N, diff));
            }

            tripletListProject.emplace_back(Eigen::Triplet<double>(i, i, -6));
            tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i, -(1+6*visc)));
            tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i, -(1+6*diff)));
        }
        laplacianProject.setFromTriplets(tripletListProject.begin(), tripletListProject.end());
        laplacianViscosity.setFromTriplets(tripletListViscosity.begin(), tripletListViscosity.end());
        laplacianDiffuse.setFromTriplets(tripletListDiffuse.begin(), tripletListDiffuse.end());
    };
};
