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
    GAUSS_SEIDEL,
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

struct Laplacian
{
    Eigen::SparseMatrix<double> A; 
    Eigen::VectorXd diag;
    Eigen::VectorXd plusi;
    Eigen::VectorXd plusj;
    Eigen::VectorXd plusk;
    Eigen::VectorXd precon;
};

struct Fluid3D
{
    Entity entity;
    double viscosity = 1.500;
    double diffusion = 0.15;
    double dt = 0.00005;
    std::uint16_t N = 128;
    Solver solver = GAUSS_SEIDEL;
    Advection advection = SEMI_LAGRANGIAN;
    bool is2D = true;

    std::vector<double> velocityFieldX = {};
    std::vector<double> velocityFieldY = {};
    std::vector<double> velocityFieldZ = {};

    std::vector<double> velocityFieldPrevX = {};
    std::vector<double> velocityFieldPrevY = {};
    std::vector<double> velocityFieldPrevZ = {};

    std::vector<double> substanceField = {};
    std::vector<double> substanceFieldPrev = {};

    std::vector<double> substanceField2 = {};
    std::vector<double> substanceField2Prev = {};

    Laplacian laplacianProject {};
    Laplacian laplacianViscosity {};
    Laplacian laplacianDiffuse {};

    std::uint32_t IX(const std::uint32_t x, const std::uint32_t y, const std::uint32_t z) const
    { 
        return is2D ? x + y * (N+2) : x + y * (N+2) + z * (N+2)*(N+2);
    };

    void setAMatricese(Laplacian& laplacian, std::uint32_t minus) const
    {
        std::uint64_t N3 = N*N*N;
        std::uint64_t N2 = N*N;
        laplacian.diag = Eigen::VectorXd::Zero(N3-minus);
        laplacian.plusi = Eigen::VectorXd::Zero(N3-minus);
        laplacian.plusj = Eigen::VectorXd::Zero(N3-minus);
        laplacian.plusk = Eigen::VectorXd::Zero(N3-minus);
        #pragma omp parallel for
        for (std::uint32_t n = 0; n < N3; ++n)
        {
            const std::uint32_t m = n % N2;
            const std::uint16_t i = m % N;
            const std::uint16_t j = m / N;
            const std::uint16_t k = n / N2;
            if (n < N3-minus)
            {
                laplacian.diag.coeffRef(n) = laplacian.A.coeff(n, n);
                if (n+1 < N3-minus && i+1 < N)
                {
                    laplacian.plusi[n] = laplacian.A.coeff(n, (i+1)+j*N+k*N2);
                }
                if (n+N < N3-minus)
                {
                    laplacian.plusj[n] = laplacian.A.coeff(n, i+(j+1)*N+k*N2);
                }
                if (n+N2 < N3-minus)
                {
                    laplacian.plusk[n] = laplacian.A.coeff(n, i+j*N+(k+1)*N2);
                }
            }
        }
    }

    void setPrecon(Laplacian& A, std::uint32_t minus) const
    {
        std::uint64_t N3 = N*N*N;
        std::uint64_t N2 = N*N;
        A.precon = Eigen::VectorXd::Zero(N3-minus);
        #pragma omp parallel for
        for (std::uint32_t n = 0; n < N3; ++n)
        {
            const std::uint32_t m = n % N2;
            const std::uint16_t i = m % N;
            const std::uint16_t j = m / N;
            const std::uint16_t k = n / N2;
            if (n < N3-minus)
            {
                std::uint32_t indmi = (i-1)+j*N+k*N2;
                std::uint32_t indmj = i+(j-1)*N+k*N2;
                std::uint32_t indmk = i+j*N+(k-1)*N2;

                double a = 0;
                double b = 0;
                double c = 0;

                double i0 = 0;
                double i1 = 0;
                double i2 = 0;
                double i3 = 0;
                double j0 = 0;
                double j1 = 0;
                double j2 = 0;
                double j3 = 0;
                double k0 = 0;
                double k1 = 0;
                double k2 = 0;
                double k3 = 0;

                if (i > 0)
                {
                    a = std::pow(A.plusi.coeff(indmi) * A.precon.coeff(indmi), 2);
                    i0 = A.plusi.coeff(indmi);
                    i1 = A.plusj.coeff(indmi);
                    i2 = A.plusk.coeff(indmi);
                    i3 = std::pow(A.precon.coeff(indmi), 2);
                }
                if (j > 0)
                {
                    b = std::pow(A.plusj.coeff(indmj) * A.precon.coeff(indmj), 2);
                    j0 = A.plusj.coeff(indmj);
                    j1 = A.plusi.coeff(indmj);
                    j2 = A.plusk.coeff(indmj);
                    j3 = std::pow(A.precon.coeff(indmj), 2);
                }
                if (k > 0)
                {
                    c = std::pow(A.plusk.coeff(indmk) * A.precon.coeff(indmk), 2);
                    k0 = A.plusk.coeff(indmk);
                    k1 = A.plusi.coeff(indmk);
                    k2 = A.plusj.coeff(indmk);
                    k3 = std::pow(A.precon.coeff(indmk), 2);
                }

                double e = A.diag.coeff(n) - a - b - c 
                    - 0.97 * (
                            i0 * (i1 + i2) * i3
                        +   j0 * (j1 + j2) * j3
                        +   k0 * (k1 + k2) * k3
                    );

                if (e < 0.25 * A.diag.coeff(n))
                {
                    e = A.diag.coeff(n);
                }
                A.precon.coeffRef(n) = 1/std::sqrt(e);
            }
        }
    }

    void init()
    {
        std::uint64_t N3 = N*N*N;
        std::uint64_t N32 = (N+2)*(N+2)*(N+2);
        velocityFieldX.reserve(N32);
        for (std::uint64_t i = 0; i < N32; ++i)
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
        substanceField2 = velocityFieldX;
        substanceField2Prev = velocityFieldX;

        double visc = dt * viscosity * N;
        double diff = dt * diffusion * N;

        const std::uint32_t minus = 1;

        Eigen::SparseMatrix<double> AProject    = Eigen::SparseMatrix<double>(N3-minus, N3-minus);
        Eigen::SparseMatrix<double> AViscosity  = Eigen::SparseMatrix<double>(N3, N3);
        Eigen::SparseMatrix<double> ADiffuse    = Eigen::SparseMatrix<double>(N3, N3);

        std::vector<Eigen::Triplet<double>> tripletListProject;
        std::vector<Eigen::Triplet<double>> tripletListViscosity;
        std::vector<Eigen::Triplet<double>> tripletListDiffuse;

        for(std::uint64_t i = 0; i < N3; ++i)
        {
            if (i+1 < N3 && static_cast<std::uint16_t>(i%N) != N-1)
            {
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i+1, -visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i+1, -diff));
            }
            if (i+1 < N3-minus && static_cast<std::uint16_t>(i%N) != N-1)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i+1, -1));
            }
            if (i+N < N3 && (i+N)%(N*N) >= N)
            {
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i+N, -visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i+N, -diff));
            }
            if (i+N < N3-minus && (i+N)%(N*N) >= N)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i+N, -1));
            }
            if (i+N*N < N3)
            {
                tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i+N*N, -visc));
                tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i+N*N, -diff));
            }
            if (i+N*N < N3-minus)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i+N*N, -1));
            }
            if (i < N3-minus)
            {
                tripletListProject.emplace_back(Eigen::Triplet<double>(i, i, 6));
            }
            tripletListViscosity.emplace_back(Eigen::Triplet<double>(i, i, 1+6*visc));
            tripletListDiffuse.emplace_back(Eigen::Triplet<double>(i, i, 1+6*diff));
        }
        AProject.setFromTriplets(tripletListProject.begin(), tripletListProject.end());
        AViscosity.setFromTriplets(tripletListViscosity.begin(), tripletListViscosity.end());
        ADiffuse.setFromTriplets(tripletListDiffuse.begin(), tripletListDiffuse.end());
        AProject = AProject+Eigen::SparseMatrix<double>(AProject.transpose());
        AViscosity = AViscosity+Eigen::SparseMatrix<double>(AViscosity.transpose());
        ADiffuse = ADiffuse+Eigen::SparseMatrix<double>(ADiffuse.transpose());
        for(std::uint32_t i = 0; i < N3; ++i)
        {
            if (i < N3-minus)
            {
                AProject.coeffRef(i, i) *= 0.5;
            }
            AViscosity.coeffRef(i, i) *= 0.5;
            ADiffuse.coeffRef(i, i) *= 0.5;
        }

        laplacianProject.A = AProject;
        setAMatricese(laplacianProject, minus);
        setPrecon(laplacianProject, minus);

        laplacianViscosity.A = AViscosity;
        setAMatricese(laplacianViscosity, 0);
        setPrecon(laplacianViscosity, 0);

        laplacianDiffuse.A = ADiffuse;
        setAMatricese(laplacianDiffuse, 0);
        setPrecon(laplacianDiffuse, 0);
    };
};
