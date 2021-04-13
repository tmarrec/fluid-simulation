#include "BasicEntities.h"
#include "Components.h"
#include <cstdint>

Mesh BasicEntities::_vector;
Mesh BasicEntities::_plane;
Mesh BasicEntities::_cube;

void BasicEntities::initBasicEntities(std::shared_ptr<Renderer> renderer)
{
    _vector =
    {
        .vertices =
        {
            0.0f, 0.0f, 0.0f,
            1.0f, 0.0f, 0.0f,
        },
        .normals =
        {
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f,
        },
        .indices =
        {
            0,  1
        },
        .renderMode = LINES,
    };
    renderer->initMesh(_vector);

    _plane =
    {
        .vertices =
        {
            -0.5f, 0.0f, -0.5f,
            0.5f, 0.0f, -0.5f,
            0.5f, 0.0f, 0.5f,
            -0.5f, 0.0f, 0.5f,
        },
        .normals =
        {
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,
        },
        .indices =
        {
            0,  1,  2,  0,  2,  3,
        },
    };

    _cube =
    {
        .vertices =
        {
            -0.5f, -0.5f, 0.5f,
            0.5f, -0.5f, 0.5f,
            0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,

            0.5f, 0.5f, 0.5f,
            0.5f, 0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, 0.5f,

            -0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, 0.5f, -0.5f,
            -0.5f, 0.5f, -0.5f,

            -0.5f, -0.5f, -0.5f,
            -0.5f, -0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, -0.5f,
                            
            0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, -0.5f,
            0.5f, 0.5f, -0.5f,

            -0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, 0.5f,
            -0.5f, -0.5f, 0.5f,
        },
        .normals =
        {
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,

            1.0f, 0.0f, 0.0f,
            1.0f, 0.0f, 0.0f,
            1.0f, 0.0f, 0.0f,
            1.0f, 0.0f, 0.0f,

            0.0f, 0.0f, -1.0f,
            0.0f, 0.0f, -1.0f,
            0.0f, 0.0f, -1.0f,
            0.0f, 0.0f, -1.0f,

            -1.0f, 0.0f, 0.0f,
            -1.0f, 0.0f, 0.0f,
            -1.0f, 0.0f, 0.0f,
            -1.0f, 0.0f, 0.0f,

            0.0f, 1.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 1.0f, 0.0f,

            0.0f, -1.0f, 0.0f,
            0.0f, -1.0f, 0.0f,
            0.0f, -1.0f, 0.0f,
            0.0f, -1.0f, 0.0f,
        },
        .indices =
        {
            0,  1,  2,  0,  2,  3,
            3,  6,  5,  4,  7,  6,
            5,  9, 11,  11, 9,  8,
            15, 12, 13, 13, 14, 15,
            17, 16, 19, 17, 19, 18,
            20, 21, 23, 21, 22, 23
        },
    };
    renderer->initMesh(_cube);
}

void BasicEntities::addTransform(Entity& entity, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
{
    gCoordinator.AddComponent(entity, Transform
    {
        .position = position,
        .rotation = rotation,
        .scale = scale
    });
}

void BasicEntities::addPlane(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
{
    auto entity = gCoordinator.CreateEntity();
    addTransform(entity, position, rotation, scale);

    gCoordinator.AddComponent(entity, _plane);

    Shader shaderProgram {};
    shaderProgram.setVert("shaders/vert.vert");
    shaderProgram.setFrag("shaders/frag.frag");

    gCoordinator.AddComponent(entity, Material
    {
        .shader = shaderProgram
    });
}

void BasicEntities::addCube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
{
    auto entity = gCoordinator.CreateEntity();
    addTransform(entity, position, rotation, scale);

    gCoordinator.AddComponent(entity, _cube);

    Shader shaderProgram {};
    shaderProgram.setVert("shaders/vert.vert");
    shaderProgram.setFrag("shaders/frag.frag");

    gCoordinator.AddComponent(entity, Material
    {
        .shader = shaderProgram
    });
}

void BasicEntities::addLineCube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
{
    auto entity = gCoordinator.CreateEntity();
    addTransform(entity, position, rotation, scale);

    gCoordinator.AddComponent(entity, Mesh
    {
        .vertices =
        {
            -0.5f, -0.5f, 0.5f,
            0.5f, -0.5f, 0.5f,
            0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,

            -0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, 0.5f, -0.5f,
            -0.5f, 0.5f, -0.5f,
        },
        .normals =
        {
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f,

            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f,
        },
        .indices =
        {
            0, 1,
            1, 2,
            2, 3,
            3, 0,

            4, 5,
            5, 6,
            6, 7,
            7, 4,

            2, 6,
            1, 5,
            0, 4,
            3, 7
        },
        .renderMode = LINES,
    });

    Shader shaderProgram {};
    shaderProgram.setVert("shaders/vert.vert");
    shaderProgram.setFrag("shaders/frag.frag");

    gCoordinator.AddComponent(entity, Material
    {
        .shader = shaderProgram
    });

}

Entity BasicEntities::addVector(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
{
    auto entity = gCoordinator.CreateEntity();

    addTransform(entity, position, rotation, scale);

    gCoordinator.AddComponent(entity, Mesh
    {
        .vertices =
        {
            0, 0, 0,
            1, 0, 0,
        },
        .normals =
        {
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f,
        },
        .indices =
        {
            0,  1
        },
        .renderMode = LINES,
    });

    Shader shaderProgram {};
    shaderProgram.setVert("shaders/vert.vert");
    shaderProgram.setFrag("shaders/frag.frag");

    gCoordinator.AddComponent(entity, Material
    {
        .shader = shaderProgram
    });
    return entity;
}

void BasicEntities::addFluid2D(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
{
    auto entity = gCoordinator.CreateEntity();
    addTransform(entity, position, rotation, scale);

    /*
    gCoordinator.AddComponent(entity, Mesh
    {
        .vertices =
        {
            -0.5f, 0.0f, -0.5f,
            0.5f, 0.0f, -0.5f,
            0.5f, 0.0f, 0.5f,
            -0.5f, 0.0f, 0.5f,
        },
        .normals =
        {
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,
        },
        .indices =
        {
            0,  1,  2,  0,  2,  3,
        },
        .renderMode = LINES
    });
    */

    Shader shaderProgram {};
    shaderProgram.setVert("shaders/vert.vert");
    shaderProgram.setFrag("shaders/frag.frag");

    gCoordinator.AddComponent(entity, Material
    {
        .shader = shaderProgram
    });

    Fluid2D fluid =
    {
        .entities = {},
        .velocityField = {},
        .velocityFieldPrev = {},
        .viscosity = 1.0f,
        .dt = 0.001f,
    };

    for (std::uint32_t i = 0; i < fluid.N+2; ++i)
    {
        for (std::uint32_t j = 0; j < fluid.N+2; ++j)
        {
            fluid.velocityField.emplace_back(glm::vec3{0,0,0});
            fluid.velocityFieldPrev.emplace_back(glm::vec3{0,0,0});

            glm::vec3 cellpos = position+(static_cast<float>(i)*glm::vec3{1.0f, 0, 0})+(static_cast<float>(j)*glm::vec3{0, 0, 1.0f})-((static_cast<float>(fluid.N+2)/2)*glm::vec3{1.0f, 0, 1.0f})+glm::vec3{0.5f, 0, 0.5f};
            fluid.entities.emplace_back(addVector(cellpos, {0,0,0}, {0.5f,0.5f,0.5f}));
        }
    }

    gCoordinator.AddComponent(entity, fluid);
}
