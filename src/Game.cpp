#include "Game.h"
#include "BasicEntities.h"

Coordinator gCoordinator;

void Game::initECS()
{
    gCoordinator.Init();
    _entities.reserve(MAX_ENTITIES);

    gCoordinator.RegisterComponent<Transform>();
    gCoordinator.RegisterComponent<Mesh>();
    gCoordinator.RegisterComponent<Camera>();
    gCoordinator.RegisterComponent<Material>();
    gCoordinator.RegisterComponent<Fluid3D>();

    _meshRendererSys = gCoordinator.RegisterSystem<MeshRenderer>();
    Camera camera
    {
        .yaw = 227,
        .pitch = -37,
        .speed = 0.1f,
        .FOV = 60,
        .transform = Transform
            {
                .position = {6.05, 5.84, 6.28},
                .rotation = {0, 0, 0},
                .scale = {1, 1, 1}
            }
    };
    _meshRendererSys->init(std::make_shared<Renderer>(_renderer), camera);
    _fluidsSys = gCoordinator.RegisterSystem<Fluids>();
    _fluidsSys->init(std::make_shared<Renderer>(_renderer));

    Signature signatureMeshRenderer;
    signatureMeshRenderer.set(gCoordinator.GetComponentType<Transform>());
    signatureMeshRenderer.set(gCoordinator.GetComponentType<Mesh>());
    signatureMeshRenderer.set(gCoordinator.GetComponentType<Material>());
    gCoordinator.SetSystemSignature<MeshRenderer>(signatureMeshRenderer);

    Signature signatureFluids;
    signatureFluids.set(gCoordinator.GetComponentType<Transform>());
    signatureFluids.set(gCoordinator.GetComponentType<Fluid3D>());
    gCoordinator.SetSystemSignature<Fluids>(signatureFluids);
}

void Game::run(WindowInfos windowInfos)
{
	_window->init(windowInfos);
	_renderer.init(_window);
    initECS();
    BasicEntities::initBasicEntities(std::make_shared<Renderer>(_renderer));

    BasicEntities::addFluid3D(glm::vec3{0,0,0}, glm::vec3{0,0,0}, glm::vec3{10,10,10});

    //BasicEntities::addCube(glm::vec3{5,0,0}, glm::vec3{0,0,0}, glm::vec3{5,5,5});

	mainLoop();
}

void Game::mainLoop()
{
    [[maybe_unused]] float dt = 0.0f;
#ifdef DEBUG_GUI
    std::deque<float> dtMean;
    std::deque<float> fluidTime;
    std::deque<float> renderTime;
    std::deque<float> inputTime;
    _fluidsSys->fluidSetupDebug();
#endif
	while (!_window->windowShouldClose())
	{
        auto startTime = std::chrono::high_resolution_clock::now();

#ifdef DEBUG_GUI
        auto fluidTimeStart = std::chrono::high_resolution_clock::now();
#endif
        _fluidsSys->update();
#ifdef DEBUG_GUI
        auto fluidTimeStop = std::chrono::high_resolution_clock::now();
#endif

#ifdef DEBUG_GUI
        _renderer.beginImgui();
        _fluidsSys->fluidDebugTool();
        float meanFluidTime = 0.0f;
        for (const auto& t : fluidTime)
            meanFluidTime += t;
        meanFluidTime /= TIME_ECHANT_NB;
        float meanRenderTime = 0.0f;
        for (const auto& t : renderTime)
            meanRenderTime += t;
        meanRenderTime /= TIME_ECHANT_NB;
        float meanInputTime = 0.0f;
        for (const auto& t : inputTime)
            meanInputTime += t;
        meanInputTime /= TIME_ECHANT_NB;
        float meanDtTime = 0.0f;
        for (const auto& t : dtMean)
            meanDtTime += t;
        meanDtTime /= TIME_ECHANT_NB;
        _renderer.debugGUI(meanDtTime, meanFluidTime, meanRenderTime, meanInputTime);
        auto renderTimeStart = std::chrono::high_resolution_clock::now();
#endif
        _renderer.prePass();
        _meshRendererSys->update();
#ifdef DEBUG_GUI
        _renderer.endImgui();
#endif
        _renderer.endPass();
        _window->swapBuffers();

#ifdef DEBUG_GUI
        auto renderTimeStop = std::chrono::high_resolution_clock::now();
        auto inputTimeStart = std::chrono::high_resolution_clock::now();
#endif
		_window->pollEvents();
#ifdef DEBUG_GUI
        auto inputTimeStop = std::chrono::high_resolution_clock::now();
#endif
        auto stopTime = std::chrono::high_resolution_clock::now();
        dt = std::chrono::duration<float, std::chrono::seconds::period>(stopTime - startTime).count();
        //std::cout << 1.0f/dt << std::endl;

#ifdef DEBUG_GUI
        fluidTime.emplace_back(std::chrono::duration<float, std::chrono::seconds::period>(fluidTimeStop - fluidTimeStart).count());
        renderTime.emplace_back(std::chrono::duration<float, std::chrono::seconds::period>(renderTimeStop - renderTimeStart).count());
        inputTime.emplace_back(std::chrono::duration<float, std::chrono::seconds::period>(inputTimeStop - inputTimeStart).count());
        dtMean.emplace_back(dt);
        if (fluidTime.size() > TIME_ECHANT_NB)
            fluidTime.pop_front();
        if (renderTime.size() > TIME_ECHANT_NB)
            renderTime.pop_front();
        if (inputTime.size() > TIME_ECHANT_NB)
            inputTime.pop_front();
        if (dtMean.size() > TIME_ECHANT_NB)
            dtMean.pop_front();
#endif
	}
}

