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
    float absorption = 100.0f;
    glm::vec3 lightIntensity = glm::vec3(1.0, 1.0, 1.0);
    std::uint32_t texture = 0;
    std::uint32_t TBO = 0;
    std::vector<float> texCoords = {};
};

struct Fluid3D
{
    Entity entity;
    double viscosity = 0.2;
    double diffusion = 0.0;
    double dt = 0.0001;
    std::uint32_t N = 24;

    std::vector<double> velocityFieldX = {};
    std::vector<double> velocityFieldY = {};
    std::vector<double> velocityFieldZ = {};

    std::vector<double> velocityFieldPrevX = {};
    std::vector<double> velocityFieldPrevY = {};
    std::vector<double> velocityFieldPrevZ = {};

    std::vector<double> substanceField = {};
    std::vector<double> substanceFieldPrev = {};

    Eigen::SparseMatrix<double> laplacianProject {};
    Eigen::SparseMatrix<double> laplacianViscosity {};
    Eigen::SparseMatrix<double> laplacianDiffuse {};

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

        std::uint32_t n = N*N*N;

        laplacianProject = Eigen::SparseMatrix<double>(n-1, n-1);
        laplacianViscosity = Eigen::SparseMatrix<double>(n, n);
        laplacianDiffuse = Eigen::SparseMatrix<double>(n, n);

        laplacianProject.reserve(Eigen::VectorXi::Constant(n-1, 7));
        laplacianViscosity.reserve(Eigen::VectorXi::Constant(n, 7));
        laplacianDiffuse.reserve(Eigen::VectorXi::Constant(n, 7));

        std::vector<Eigen::Triplet<double>> tripletListProject;
        std::vector<Eigen::Triplet<double>> tripletListViscosity;
        std::vector<Eigen::Triplet<double>> tripletListDiffuse;

        tripletListProject.reserve(7*(n-1));
        tripletListViscosity.reserve(7*(n));
        tripletListDiffuse.reserve(7*(n));

        std::uint32_t entities = 0;
        for(std::uint32_t i = 0; i < n; ++i)
        {
            if (i > 0 && (i-1)%N != 0)
            {
                if (i < n-1)
                {
                    tripletListProject.emplace_back(Eigen::Triplet<double>(i, i-1, 1));
                }
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i-1, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i-1, diff));
                entities++;
            }
            if (i >= N)
            {
                if (i < n-1)
                {
                    tripletListProject.emplace_back(Eigen::Triplet<double>(i, i-N, 1));
                }
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i-N, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i-N, diff));
                entities++;
            }
            if (i >= N * N)
            {
                if (i < n-1)
                {
                    tripletListProject.emplace_back(Eigen::Triplet<double>(i, i-N*N, 1));
                }
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i-N*N, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i-N*N, diff));
                entities++;
            }
            if (i < N*N*N-1 && (i+1)%N != 0)
            {
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i+1, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i+1, diff));
            }
            if (i < N*N*N-1-1 && (i+1)%N != 0)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i+1, 1));
                entities++;
            }
            if (i + N < N*N*N-1)
            {
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i+N, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i+N, diff));
            }
            if (i + N < N*N*N-1-1)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i+N, 1));
                entities++;
            }
            if (i + N*N < N*N*N-1)
            {
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i+N*N, visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i+N*N, diff));
            }
            if (i + N*N < N*N*N-1-1)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i+N*N, 1));
                entities++;
            }
            if (i < n-1)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i, -6));
            }
            tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i, -(1+6*visc)));
            tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i, -(1+6*diff)));
        }
        laplacianProject.setFromTriplets(tripletListProject.begin(), tripletListProject.end());
        laplacianViscosity.setFromTriplets(tripletListViscosity.begin(), tripletListViscosity.end());
        laplacianDiffuse.setFromTriplets(tripletListDiffuse.begin(), tripletListDiffuse.end());

    };
};
