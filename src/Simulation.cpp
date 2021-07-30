#include "Simulation.h"

void Simulation::run(WindowInfos windowInfos)
{
    //_renderer = new Renderer(windowInfos);
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
        if (it == 64)
            exit(0);
            */
        //std::cout << "== Iteration " << iterations << " ==" << std::endl;
            
        auto startTime = std::chrono::high_resolution_clock::now();

        _fluid.update(it);

	    _renderer.initTexture2D(_fluid.texture(), _fluidRenderer.material.texture);
        updateMeshVec();
        updateMeshGrid();

        _renderer.prePass();


        _renderer.applyMaterial(_fluidRenderer.material, _camera, _fluidRenderer.transform);
        _renderer.drawMesh(_fluidRenderer.mesh);

        _renderer.applyMaterial(_fluidRenderer.materialVec, _camera, _fluidRenderer.transform);
        _renderer.drawMesh(_fluidRenderer.meshVec);

        _renderer.applyMaterial(_fluidRenderer.materialGrid, _camera, _fluidRenderer.transform);
        _renderer.drawMesh(_fluidRenderer.meshGrid);

        _renderer.endPass();

        //_renderer.writeImg(it);
        _window.swapBuffers();

		_window.pollEvents();
        auto stopTime = std::chrono::high_resolution_clock::now();
        dt = std::chrono::duration<float, std::chrono::seconds::period>(stopTime - startTime).count();
        //std::cout << "Done in " << dt << " sec" << std::endl << std::endl;

        it++;
	}
}

void Simulation::updateMeshGrid()
{
    Mesh& mesh = _fluidRenderer.meshGrid;
    const std::uint16_t& N = _fluid.N();

    float z = 0.001f;
    std::uint64_t it = 0;

    mesh.vertices.clear();
    mesh.indices.clear();
    for (float i = 0; i <= N; ++i)
    {
        // Arrow line drawing
        glm::vec2 A = { (i/N)-0.5f, -0.5f };
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

void Simulation::updateMeshVec()
{
    Mesh& mesh = _fluidRenderer.meshVec;
    const std::uint16_t& N = _fluid.N();
    const std::vector<double>& X = _fluid.X();
    const std::vector<double>& Y = _fluid.Y();

    float z = 0.001f;
    std::uint64_t it = 0;
    float reduce = 600.0f;

    mesh.vertices.clear();
    mesh.indices.clear();
    for (float j = 0; j < N; ++j)
    {
        for (float i = 0; i < N+1; ++i)
        {
            if (_fluid.isCellActive(i,j))
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
            }
        }
    }

    for (float j = 0; j < N+1; ++j)
    {
        for (float i = 0; i < N; ++i)
        {
            if (_fluid.isCellActive(i,j))
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
            }
        }
    }
    /*
    for (float i = 0; i < N+1; ++i)
    {
        for (float j = 0; j < N+1; ++j)
        {
            // Arrow line drawing
            glm::vec2 A = { ((i+0.5)/N)-0.5, ((j+0.5)/N)-0.5 };
            mesh.vertices.emplace_back(A.x);
            mesh.vertices.emplace_back(z);
            mesh.vertices.emplace_back(A.y);

            float u = 0.5*(X[i+j*(N+1)]+X[(i+1)+j*(N+1)]);
            float v = 0.5*(Y[i+j*(N)]+Y[i+(j+1)*(N)]);

            glm::vec2 uv = {u,v};
            uv = glm::normalize(uv);
            u = uv.x*0.5;
            v = uv.y*0.5;

            glm::vec2 B = { ((i+0.5+u)/N)-0.5, ((j+0.5+v)/N)-0.5 };

            mesh.vertices.emplace_back(B.x);
            mesh.vertices.emplace_back(z);
            mesh.vertices.emplace_back(B.y);

            // Arrow head drawing
            float size = (1.0/N)/5;
            float h = size*sqrtf(3), w = size;
            glm::vec2 U = glm::normalize(B-A);
            glm::vec2 V = glm::vec2(-U.y, U.x);
            glm::vec2 v1 = B - h*U + w*V;
            glm::vec2 v2 = B - h*U - w*V;

            mesh.vertices.emplace_back(v1.x);
            mesh.vertices.emplace_back(z);
            mesh.vertices.emplace_back(v1.y);

            mesh.vertices.emplace_back(v2.x);
            mesh.vertices.emplace_back(z);
            mesh.vertices.emplace_back(v2.y);

            mesh.indices.emplace_back(it);
            mesh.indices.emplace_back(it+1);
            mesh.indices.emplace_back(it+1);
            mesh.indices.emplace_back(it+2);
            mesh.indices.emplace_back(it+1);
            mesh.indices.emplace_back(it+3);
            it += 4;
        }
    }
    */
    _renderer.initMesh(mesh);
}

void Simulation::initSimulation()
{
    Shader shaderProgram {};
    shaderProgram.setVert("shaders/vert2D.vert");
    shaderProgram.setFrag("shaders/fluid2D.frag");
    Material material = 
    {
        .shader = shaderProgram,
        .is2D = true,
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
        .is2D = true,
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
    _renderer.initMesh(_fluidRenderer.meshVec);

    Shader shaderProgramGrid {};
    shaderProgramGrid.setVert("shaders/vert.vert");
    shaderProgramGrid.setFrag("shaders/grid.frag");
    Material materialGrid = 
    {
        .shader = shaderProgramGrid,
        .is2D = true,
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
    _renderer.initMesh(_fluidRenderer.meshGrid);

    // Set camera front
    glm::vec3 dir;
    dir.x = cos(glm::radians(_camera.yaw))*cos(glm::radians(_camera.pitch));
    dir.y = sin(glm::radians(_camera.pitch));
    dir.z = sin(glm::radians(_camera.yaw))*cos(glm::radians(_camera.pitch));
    dir = glm::normalize(dir);
    _camera.front = dir;
}
