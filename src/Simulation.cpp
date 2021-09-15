#include "Simulation.h"
#include "config.h"
#include <glm/gtx/string_cast.hpp>

void Simulation::run(WindowInfos windowInfos)
{
	_window.init(windowInfos);
	_renderer.init(windowInfos);
    initSimulation();
	mainLoop();
}

void Simulation::mainLoop()
{
    [[maybe_unused]] float dt = 0.0f;
    std::uint64_t it = 0;

	while (!_window.windowShouldClose())
	{
        /*
        std::cout << _camera.yaw << std::endl;
        std::cout << _camera.pitch << std::endl;
        std::cout << glm::to_string(_camera.transform.position) << std::endl;
        */
        if (it == 512)
        {
            exit(0);
        }
        //std::cout << "== Iteration " << iterations << " ==" << std::endl;
            
        auto startTime = std::chrono::high_resolution_clock::now();

        _fluid.update(it);

        handleInputs();
        setCameraDir();
        if (Config::dim == 2)
        {
	        _renderer.initTexture2D(_fluid.texture(), _fluidRenderer.material.texture);
            updateMeshVec();
            updateMeshGrid();
            updateMeshGridBorder();
        }
        else if (Config::dim == 3)
        {
	        _renderer.initTexture3D(_fluid.texture(), _fluidRenderer.material.texture);
        }

        _renderer.prePass();

        _renderer.applyMaterial(_fluidRenderer.material, _camera, _fluidRenderer.transform);
        _renderer.drawMesh(_fluidRenderer.mesh);

        if (Config::dim == 2)
        {
            /*
            _renderer.applyMaterial(_fluidRenderer.materialVec, _camera, _fluidRenderer.transform);
            _renderer.drawMesh(_fluidRenderer.meshVec);

            _renderer.applyMaterial(_fluidRenderer.materialGrid, _camera, _fluidRenderer.transform);
            _renderer.drawMesh(_fluidRenderer.meshGrid);
            */

            _renderer.setLineWidth(2);
            _renderer.applyMaterial(_fluidRenderer.materialGridBorder, _camera, _fluidRenderer.transform);
            _renderer.drawMesh(_fluidRenderer.meshGridBorder);
            _renderer.setLineWidth(1);
        }

        _renderer.endPass();
        if (Config::dim == 3)
        {
            _renderer.raymarchPass();
            marchingCube.run(_fluid.surface(), it);
        }

        if (Config::exportFrames)
        {
            _renderer.writeImg(it);
        }
        _window.swapBuffers();
        _window.pollEvents();

        auto stopTime = std::chrono::high_resolution_clock::now();
        dt = std::chrono::duration<float, std::chrono::seconds::period>(stopTime - startTime).count();
        //std::cout << "Done in " << dt << " sec" << std::endl << std::endl;

        it++;
	}

    _renderer.freeMesh(_fluidRenderer.mesh);
    _renderer.freeMesh(_fluidRenderer.meshGrid);
    _renderer.freeMesh(_fluidRenderer.meshGridBorder);
    _renderer.freeMesh(_fluidRenderer.meshVec);
}

void Simulation::updateMeshGrid()
{
    Mesh& mesh = _fluidRenderer.meshGrid;
    const std::uint16_t& N = Config::N;

    float z = 0.001f;
    std::uint64_t it = 0;

    mesh.vertices.clear();
    mesh.indices.clear();
    for (float i = 0; i <= N; ++i)
    {
        // Arrow line drawing
        glm::vec2 A = { (i/N)-0.5f, -1.5f };
        mesh.vertices.emplace_back(A.x);
        mesh.vertices.emplace_back(z);
        mesh.vertices.emplace_back(A.y);

        glm::vec2 B = A + glm::vec2{0.0f, 1.0f};
        mesh.vertices.emplace_back(B.x);
        mesh.vertices.emplace_back(z);
        mesh.vertices.emplace_back(B.y);


        mesh.vertices.emplace_back(A.y);
        mesh.vertices.emplace_back(z);
        mesh.vertices.emplace_back(A.x);

        mesh.vertices.emplace_back(B.y);
        mesh.vertices.emplace_back(z);
        mesh.vertices.emplace_back(B.x);

        mesh.indices.emplace_back(it);
        mesh.indices.emplace_back(it+1);
        mesh.indices.emplace_back(it+2);
        mesh.indices.emplace_back(it+3);
        it += 4;
    }

    _renderer.initMesh(mesh);
}

void Simulation::updateMeshGridBorder()
{
    Mesh& mesh = _fluidRenderer.meshGridBorder;
    float z = 0.0011f;

    mesh.vertices.clear();
    mesh.indices.clear();

    mesh.vertices.emplace_back(-0.5f);
    mesh.vertices.emplace_back(z);
    mesh.vertices.emplace_back(-0.5f);

    mesh.vertices.emplace_back(0.5f);
    mesh.vertices.emplace_back(z);
    mesh.vertices.emplace_back(-0.5f);

    mesh.vertices.emplace_back(0.5f);
    mesh.vertices.emplace_back(z);
    mesh.vertices.emplace_back(0.5f);

    mesh.vertices.emplace_back(-0.5f);
    mesh.vertices.emplace_back(z);
    mesh.vertices.emplace_back(0.5f);

    mesh.indices.emplace_back(0);
    mesh.indices.emplace_back(1);

    mesh.indices.emplace_back(1);
    mesh.indices.emplace_back(2);

    mesh.indices.emplace_back(2);
    mesh.indices.emplace_back(3);

    mesh.indices.emplace_back(3);
    mesh.indices.emplace_back(0);

    _renderer.initMesh(mesh);
}

void Simulation::updateMeshVec()
{
    Mesh& mesh = _fluidRenderer.meshVec;
    const std::uint16_t& N = Config::N;
    const std::vector<double>& X = _fluid.X();
    const std::vector<double>& Y = _fluid.Y();

    float z = 0.001f;
    std::uint64_t it = 0;
    float reduce = 18000.0f;

    mesh.vertices.clear();
    mesh.indices.clear();
    for (float j = 0; j < N; ++j)
    {
        for (float i = 0; i < N+1; ++i)
        {
            //if (_fluid.isCellActive(i,j,0))
            {
                float size = float(X[i+j*(N+1)]/reduce);

                glm::vec2 A = { (i/N)-0.5f, ((j+0.5f)/N)-0.5f };
                mesh.vertices.emplace_back(A.x);
                mesh.vertices.emplace_back(z);
                mesh.vertices.emplace_back(A.y);

                glm::vec2 B = A + glm::vec2{size, 0.0f };
                mesh.vertices.emplace_back(B.x);
                mesh.vertices.emplace_back(z);
                mesh.vertices.emplace_back(B.y);

                mesh.indices.emplace_back(it);
                mesh.indices.emplace_back(it+1);
                it += 2;

                size = float(X[(i+1)+j*(N+1)]/reduce);

                A = { ((i+1)/N)-0.5f, ((j+0.5f)/N)-0.5f };
                mesh.vertices.emplace_back(A.x);
                mesh.vertices.emplace_back(z);
                mesh.vertices.emplace_back(A.y);

                B = A + glm::vec2{size, 0.0f };
                mesh.vertices.emplace_back(B.x);
                mesh.vertices.emplace_back(z);
                mesh.vertices.emplace_back(B.y);

                mesh.indices.emplace_back(it);
                mesh.indices.emplace_back(it+1);
                it += 2;
            }
        }
    }

    for (float j = 0; j < N+1; ++j)
    {
        for (float i = 0; i < N; ++i)
        {
            //if (_fluid.isCellActive(i,j,0))
            {
                float size = float(Y[i+j*N]/reduce);

                glm::vec2 A = { ((i+0.5f)/N)-0.5f, (j/N)-0.5f };
                mesh.vertices.emplace_back(A.x);
                mesh.vertices.emplace_back(z);
                mesh.vertices.emplace_back(A.y);

                glm::vec2 B = A + glm::vec2{0.0f, size};
                mesh.vertices.emplace_back(B.x);
                mesh.vertices.emplace_back(z);
                mesh.vertices.emplace_back(B.y);

                mesh.indices.emplace_back(it);
                mesh.indices.emplace_back(it+1);
                it += 2;

                size = float(Y[i+(j+1)*N]/reduce);

                A = { ((i+0.5f)/N)-0.5f, ((j+1)/N)-0.5f };
                mesh.vertices.emplace_back(A.x);
                mesh.vertices.emplace_back(z);
                mesh.vertices.emplace_back(A.y);

                B = A + glm::vec2{0.0f, size};
                mesh.vertices.emplace_back(B.x);
                mesh.vertices.emplace_back(z);
                mesh.vertices.emplace_back(B.y);

                mesh.indices.emplace_back(it);
                mesh.indices.emplace_back(it+1);
                it += 2;
            }
        }
    }
    _renderer.initMesh(mesh);
}

void Simulation::initSimulation()
{
    Shader shaderProgram {};

    _fluidRenderer.mesh.dim = Config::dim;
    if (Config::dim == 3)
    {
        shaderProgram.setVert("shaders/vert.vert");
        shaderProgram.setFrag("shaders/fluid3D.frag");
    }
    else if (Config::dim == 2)
    {
        shaderProgram.setVert("shaders/vert2D.vert");
        shaderProgram.setFrag("shaders/fluid2D.frag");
    }
    Material material = 
    {
        .shader = shaderProgram,
        .dim = Config::dim,
        .hasTexture = true,
        .texCoords =
        {
            1.0f, 1.0f,
            1.0f, 0.0f,
            0.0f, 0.0f, 
            0.0f, 1.0f
        },
    };
    _fluidRenderer.material = material;
    _renderer.initMaterial(material);
    _renderer.initMesh(_fluidRenderer.mesh);

    Shader shaderProgramVec {};
    shaderProgramVec.setVert("shaders/vert.vert");
    shaderProgramVec.setFrag("shaders/vec.frag");
    Material materialVec = 
    {
        .shader = shaderProgramVec,
        .texCoords =
        {
            1.0f, 1.0f,
            1.0f, 0.0f,
            0.0f, 0.0f, 
            0.0f, 1.0f
        },
    };
    _fluidRenderer.materialVec = materialVec;
    _renderer.initMaterial(materialVec);

    Shader shaderProgramGrid {};
    shaderProgramGrid.setVert("shaders/vert.vert");
    shaderProgramGrid.setFrag("shaders/grid.frag");
    Material materialGrid = 
    {
        .shader = shaderProgramGrid,
        .texCoords =
        {
            1.0f, 1.0f,
            1.0f, 0.0f,
            0.0f, 0.0f, 
            0.0f, 1.0f
        },
    };
    _fluidRenderer.materialGrid = materialGrid;
    _renderer.initMaterial(materialGrid);

    Shader shaderProgramGridBorder {};
    shaderProgramGridBorder.setVert("shaders/vert.vert");
    shaderProgramGridBorder.setFrag("shaders/gridBorder.frag");
    Material materialGridBorder = 
    {
        .shader = shaderProgramGridBorder,
        .texCoords =
        {
            1.0f, 1.0f,
            1.0f, 0.0f,
            0.0f, 0.0f, 
            0.0f, 1.0f
        },
    };
    _fluidRenderer.materialGridBorder = materialGridBorder;
    _renderer.initMaterial(materialGridBorder);

    if (Config::dim == 3)
    {
        _camera = 
        {
            .yaw = -572,
            .pitch = -31,
            .speed = 1,
            .transform = Transform
            {
                .position = {14.85, 9.39, -8.22},
                .rotation = {90, 0, 0},
                .scale = {1, 1, 1}
            },
        };
    }
    else if (Config::dim == 2)
    {
        _camera = 
        {
            .yaw = -90-90-90,
            .pitch = -90,
            .speed = 1,
            .transform = Transform
            {
                .position = {0.0, 1.0, 0.0},
                .rotation = {0, 0, 0},
                .scale = {1, 1, 1}
            },
        };
    }
}

void Simulation::setCameraDir()
{
    // Set camera front
    glm::vec3 dir;
    dir.x = cos(glm::radians(_camera.yaw))*cos(glm::radians(_camera.pitch));
    dir.y = sin(glm::radians(_camera.pitch));
    dir.z = sin(glm::radians(_camera.yaw))*cos(glm::radians(_camera.pitch));
    dir = glm::normalize(dir);
    _camera.front = dir;
}

void Simulation::handleInputs()
{
    if (Input::keyIsDown(GLFW_KEY_W))
    {
        _camera.transform.position += _camera.front * _camera.speed;
    }
    if (Input::keyIsDown(GLFW_KEY_A))
    {
        _camera.transform.position -= glm::normalize(glm::cross(_camera.front, _camera.up)) * _camera.speed;
    }
    if (Input::keyIsDown(GLFW_KEY_S))
    {
        _camera.transform.position -= _camera.front * _camera.speed;
    }
    if (Input::keyIsDown(GLFW_KEY_D))
    {
        _camera.transform.position += glm::normalize(glm::cross(_camera.front, _camera.up)) * _camera.speed;
    }
    if (Input::keyIsDown(GLFW_KEY_LEFT_SHIFT))
    {
        _camera.transform.position.y += _camera.speed;
    }
    if (Input::keyIsDown(GLFW_KEY_LEFT_CONTROL))
    {
        _camera.transform.position.y -= _camera.speed;
    }

    float sensitivity = 1.0f;
    _camera.yaw += Input::mouseOffsetX*sensitivity;
    _camera.pitch = glm::clamp(_camera.pitch - Input::mouseOffsetY*sensitivity, -89.9f, 89.9f);

    glm::vec3 dir;
    dir.x = cos(glm::radians(_camera.yaw))*cos(glm::radians(_camera.pitch));
    dir.y = sin(glm::radians(_camera.pitch));
    dir.z = sin(glm::radians(_camera.yaw))*cos(glm::radians(_camera.pitch));
    dir = glm::normalize(dir);
    _camera.front = dir;

    Input::updateMouseMovements();
}
