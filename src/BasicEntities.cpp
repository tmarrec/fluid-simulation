#include "BasicEntities.h"

Mesh BasicEntities::_vector;
Mesh BasicEntities::_plane;
Mesh BasicEntities::_cube;
std::shared_ptr<Renderer> BasicEntities::_renderer;

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

    Mesh mesh
    {
        .vertices =
        {
            -0.5f, 0.1f, -0.5f,
            0.5f, 0.1f, -0.5f,
            0.5f, 0.1f, 0.5f,
            -0.5f, 0.1f, 0.5f,
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
    gCoordinator.AddComponent(entity, mesh);

    Shader shaderProgram {};
    shaderProgram.setVert("shaders/vert.vert");
    shaderProgram.setFrag("shaders/fluid2D.frag");

    Material material = 
    {
        .shader = shaderProgram,
        .texCoords =
        {
            1.0f, 1.0f,
            1.0f, 0.0f,
            0.0f, 0.0f, 
            0.0f, 1.0f
        }
    };
    BasicEntities::_renderer->initMaterial(material);
    gCoordinator.AddComponent(entity, material);


    Fluid2D fluid =
    {
        .viscosity = 0.02f,
        .dt = 0.0001f,
        .entity = entity,
    };

    fluid.velocityFieldX.reserve((fluid.N+2)*(fluid.N+2));
    for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2); ++i)
    {
        fluid.velocityFieldX.emplace_back(0);
    }
    fluid.velocityFieldPrevX = fluid.velocityFieldX;
    fluid.velocityFieldY = fluid.velocityFieldX;
    fluid.velocityFieldPrevY = fluid.velocityFieldX;
    fluid.substanceField = fluid.velocityFieldX;
    fluid.substanceFieldPrev = fluid.velocityFieldX;

    gCoordinator.AddComponent(entity, fluid);
}

void BasicEntities::addFluid3D(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
{
    auto entity = gCoordinator.CreateEntity();
    addTransform(entity, position, rotation, scale);

    Mesh mesh
    {
        .vertices =
        {
            -0.5f, 0.1f, -0.5f,
            0.5f, 0.1f, -0.5f,
            0.5f, 0.1f, 0.5f,
            -0.5f, 0.1f, 0.5f,
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
    gCoordinator.AddComponent(entity, mesh);

    Material material = 
    {
        .noShader = true,
        .texCoords =
        {
            1.0f, 1.0f,
            1.0f, 0.0f,
            0.0f, 0.0f, 
            0.0f, 1.0f
        },
    };
    BasicEntities::_renderer->initMaterial(material);
    gCoordinator.AddComponent(entity, material);

    Fluid3D fluid =
    {
        .entity = entity,
    };
    fluid.init();

    gCoordinator.AddComponent(entity, fluid);
}
