#include "GlWidget.h"
#include "src/ui/MainWindow.h"

#include <GL/gl.h>
#include <cstdint>
#include <qevent.h>
#include <unistd.h>
#include <thread>

#include "InputManager.h"
#include "shapes.h"
#include "Renderer.h"
#include "TransformComponent.h"
#include "CameraComponent.h"
#include "DrawableComponent.h"
#include "LightComponent.h"
#include "SubdivideComponent.h"

#include "BSpline.h"
#include "BSplineTensor.h"

#include "Scene.h"

void GLAPIENTRY MessageCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam);

GlWidget::GlWidget(QWidget *parent)
: QOpenGLWidget{parent}
, _renderer { std::shared_ptr<Renderer>(new Renderer) }
, _ECS_manager { std::shared_ptr<ECS_Manager>(new ECS_Manager) }
, _InputManager { std::shared_ptr<InputManager>(new InputManager(static_cast<GlWidget*>(this))) }
{
	setFocus();
	//_currentTime = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

Entity* GlWidget::getActiveCamera() const { return _activeCamera; }	

void GlWidget::_init()
{

	_activeCamera = &_ECS_manager->addEntity();
	_activeCamera->addComponent<CameraComponent>(0.0f, 0.0f, 0.05f, 90.0f);
	_activeCamera->addComponent<TransformComponent>(glm::vec3{-5.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
	_renderer->setActiveCamera(&_activeCamera->getComponent<CameraComponent>());

	auto shader = std::make_shared<Shader>(Shader{"shaders/vert.vert", "shaders/frag.frag"});

	auto s = Scene(_renderer, _ECS_manager, "TEST.gltf");

	Cube c;

	auto & light(_ECS_manager->addEntity());
	light.addComponent<TransformComponent>(glm::vec3{0.0f, 100.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{20.0f, 20.0f, 20.0f});
	light.addComponent<DrawableComponent>(_renderer, shader, c.vertices, c.normals, c.indices, GL_TRIANGLES);
	light.addComponent<LightComponent>(_renderer, glm::vec3{1.0f, 0.0f, 0.0f}, 1.f);

	auto & light2(_ECS_manager->addEntity());
	light2.addComponent<TransformComponent>(glm::vec3{0.0f, 200.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{20.0f, 20.0f, 20.0f});
	light2.addComponent<DrawableComponent>(_renderer, shader, c.vertices, c.normals, c.indices, GL_TRIANGLES);
	light2.addComponent<LightComponent>(_renderer, glm::vec3{1.0f, 1.0f, 1.0f}, 4.f);

	/*
	auto & cube(_ECS_manager->addEntity());
	cube.addComponent<TransformComponent>(glm::vec3{-50.0f, -20.0f, -50.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	cube.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", c.vertices, c.normals, c.indices, GL_TRIANGLES);

	auto & cube2(_ECS_manager->addEntity());
	cube2.addComponent<TransformComponent>(glm::vec3{20.0f, 10.0f, 80.0f}, glm::vec3{30.0f, 10.0f, 60.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	//cube2.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", c.vertices, c.normals, c.indices, GL_TRIANGLES);

	auto & bsplineline(_ECS_manager->addEntity());
	std::vector<glm::vec3> controls = 
	{
			{0, 0, 0},
			{1, 1.5f, 0},
			{2, -1.5f, 0},
			{3, 0, 0},
			{4, 1.5f, 0},
	};
	auto bsp = BSplineShape(3, controls, true, 0.1f).shape();
	bsplineline.addComponent<TransformComponent>(glm::vec3{0.0f, 50.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	//bsplineline.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", bsp.vertices, bsp.normals, bsp.indices, GL_LINE_STRIP);

	auto & bsplineline2(_ECS_manager->addEntity());
	std::vector<glm::vec3> controls2 =
	{
			{0, 0, 0},
			{0, 2, 0},
			{1, 2, 0},
			{4, 0, 0},
			{5, 0, 0},
			{5, 2, 0},
	};
	auto bsp2 = BSplineShape(3, controls2, true, 0.01f, {0,0,0,0.25f,0.5f,0.75f,1,1,1}).shape();
	bsplineline2.addComponent<TransformComponent>(glm::vec3{-200.0f, 0.0f, 400.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	//bsplineline2.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", bsp2.vertices, bsp2.normals, bsp2.indices, GL_LINE_STRIP);



	auto & bsplinetensor(_ECS_manager->addEntity());
	std::vector<std::vector<glm::vec3>> test =
	{
		{
			{0, 0, 0},
			{1, 1.5f, 0},
			{2, -1.5f, 0},
			{3, 0, 0},
			{4, 1, 0},
		},
		{
			{0, 0, 1},
			{1, 1.5f, 1},
			{2, -1.5f, 1},
			{3, 0, 1},
			{4, 1, 1},
		},
		{
			{0, 0, 2},
			{1, 1.5f, 2},
			{2, -1.5f, 2},
			{3, 0, 2},
			{4, 1, 2},
		}
	};	
	auto bspt = BSplineTensorShape(3, test, 0.15f).shape();
	bsplinetensor.addComponent<TransformComponent>(glm::vec3{0.0f, 50.0f, -150.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	//bsplinetensor.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", bspt.vertices, bspt.normals, bspt.indices, GL_TRIANGLES);

	auto & bsplinetensorknots(_ECS_manager->addEntity());
	std::vector<std::vector<glm::vec3>> test2 =
	{
		{
			{0, 0, 0},
			{0, 2, 0},
			{1, 2, 0},
			{4, 0, 0},
			{5, 0, 0},
			{5, 2, 0},
		},
		{
			{0, 0, 2},
			{0, 2, 2},
			{1, 2, 2},
			{4, 0, 2},
			{5, 0, 2},
			{5, 2, 2},
		},
		{
			{0, 0, 4},
			{0, 2, 4},
			{1, 2, 4},
			{4, 0, 4},
			{5, 0, 4},
			{5, 2, 4},
		}
	};	
	std::vector<std::vector<float>> knots =
	{
		{0,0,0,0.25f,0.5f,0.75f,1,1,1},
		{0,0,0,0,0,0,1.0f,1,1},
		{0,0,0,0.25f,0.5f,0.75f,1,1,1},
	};
	auto bsptknots = BSplineTensorShape(3, test2, 0.004f, knots).shape();
	bsplinetensorknots.addComponent<TransformComponent>(glm::vec3{-200.0f, 0.0f, 150.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	//bsplinetensorknots.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", bsptknots.vertices, bsptknots.normals, bsptknots.indices, GL_TRIANGLES);

	auto & bsplinetensor2(_ECS_manager->addEntity());
	std::vector<std::vector<glm::vec3>> randomControls;
	for (int i = 0; i < 40; ++i) {
		std::vector<glm::vec3> temp;
		for (int j = 0; j < 40; ++j) {
			temp.push_back({j, rand()%4, i});
		}
		randomControls.emplace_back(temp);
	}
	auto bspt2 = BSplineTensorShape(4, randomControls, 0.20f).shape();
	bsplinetensor2.addComponent<TransformComponent>(glm::vec3{-1000.0f, -200.0f, -1000.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	//bsplinetensor2.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", bspt2.vertices, bspt2.normals, bspt2.indices, GL_TRIANGLES);
	

	auto & CUBE(_ECS_manager->addEntity());
	CUBE.addComponent<TransformComponent>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{150.0f, 150.0f, 150.0f});
	CUBE.addComponent<DrawableComponent>(_renderer, shader, c.vertices, c.normals, c.indices, GL_TRIANGLES);
	CUBE.addComponent<SubdivideComponent>();

	auto & floor(_ECS_manager->addEntity());
	floor.addComponent<TransformComponent>(glm::vec3{0.0f, -200.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1000.0f, 10.0f, 1000.0f});
	floor.addComponent<DrawableComponent>(_renderer, shader, c.vertices, c.normals, c.indices, GL_TRIANGLES);
	*/
}

void GlWidget::initializeGL()
{
	connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GlWidget::cleanup);

	if (!initializeOpenGLFunctions()) {
		exit(1);
	}

	glEnable(GL_DEBUG_OUTPUT);
	glEnable(GL_DEPTH_TEST);
	glDebugMessageCallback(MessageCallback, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	_renderer->initGl();
	_init();
	//std::thread t ((Test(_ECS_manager)));
	//t.detach();
}

void GlWidget::switchPolygonmode()
{
	makeCurrent();
	if (_polygonFillMode)
	{
		_polygonFillMode = false;
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else
	{
		_polygonFillMode = true;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	doneCurrent();
}

void GlWidget::paintGL()
{
	_renderer->startFrame();
	_ECS_manager->update(_deltaTime);
	_ECS_manager->draw();
	_renderer->endFrame();

	//#####################################
	// Compte les FPS chaque secondes
	std::uint64_t end_timer_fps = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	if (end_timer_fps-_start_timer_fps > 1000) {
		//_main_window->update_title_infos("FPS: " + std::to_string(_frame_count));
		std::cout << _frame_count << std::endl;
		_frame_count = 0;
		_start_timer_fps = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	}

	// Calcule le Time Delta
	if (_start_timer_frame != 0) {
		std::uint64_t end_timer_frame = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		_deltaTime = 1e-9*(end_timer_frame-_start_timer_frame);
	}
	_start_timer_frame = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	_frame_count++;
	//#####################################

	update();
}

void GlWidget::resizeGL(int w, int h)
{
	makeCurrent();
	_renderer->resizeGl(w, h);
	doneCurrent();
}
void GlWidget::keyPressEvent(QKeyEvent *event) { _InputManager->keyPressEvent(event); }
void GlWidget::keyReleaseEvent(QKeyEvent *event) { _InputManager->keyReleaseEvent(event); }
void GlWidget::mousePressEvent(QMouseEvent *event) { _InputManager->mousePressEvent(event); }
void GlWidget::mouseMoveEvent(QMouseEvent *event) { _InputManager->mouseMoveEvent(event); }

void GlWidget::cleanup()
{

}

GlWidget::~GlWidget()
{

}

void GLAPIENTRY
MessageCallback
(
	[[maybe_unused]] GLenum source,
	GLenum type,
	[[maybe_unused]] GLuint id,
	GLenum severity,
	[[maybe_unused]] GLsizei length,
	const GLchar* message,
	[[maybe_unused]] const void* userParam
)
{
	if (severity != 33387)
	{
		fprintf(stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
			(type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""), type, severity, message);
	}
}

