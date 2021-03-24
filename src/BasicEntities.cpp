#include "BasicEntities.h"

void BasicEntities::addTransform(Entity& entity, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
{
    gCoordinator.AddComponent(entity, Transform
    {
        .position = position,
        .rotation = rotation,
        .scale = scale
    });
}

void BasicEntities::addCube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
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
    });

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

void BasicEntities::addVector(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
{
    auto entity = gCoordinator.CreateEntity();

    addTransform(entity, position, rotation, scale);

    gCoordinator.AddComponent(entity, Mesh
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
    });

    Shader shaderProgram {};
    shaderProgram.setVert("shaders/vert.vert");
    shaderProgram.setFrag("shaders/frag.frag");

    gCoordinator.AddComponent(entity, Material
    {
        .shader = shaderProgram
    });
}
