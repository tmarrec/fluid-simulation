#include "GlWidget.h"
#include "ui/MainWindow.h"

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
#include "MetaballComponent.h"
#include "MarchingCubeComponent.h"

#include "BSpline.h"
#include "BSplineTensor.h"

#include "Scene.h"

void GLAPIENTRY MessageCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam);

GlWidget::GlWidget(QWidget *parent)
: QOpenGLWidget{parent}
, _ECS_manager { std::shared_ptr<ECS_Manager>(new ECS_Manager) }
, _renderer { std::shared_ptr<Renderer>(new Renderer(_ECS_manager)) }
, _InputManager { std::shared_ptr<InputManager>(new InputManager(static_cast<GlWidget*>(this), _renderer)) }
{
	setFocus();
}

void GlWidget::_initScene()
{
	auto camera = &_ECS_manager->addEntity();
	camera->addComponent<CameraComponent>(0.0f, 0.0f, 0.05f, 90.0f);
	camera->addComponent<TransformComponent>(glm::vec3{-5.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
	_renderer->setActiveCamera(&camera->getComponent<CameraComponent>());

	auto shader = std::make_shared<Shader>(Shader{"shaders/vert.vert", "shaders/frag.frag"});
	auto s = Scene(_renderer, _ECS_manager, "models/inside-cube-empty.gltf");
	Cube c;

	auto& light(_ECS_manager->addEntity());
	light.addComponent<TransformComponent>(glm::vec3{-8.0f, 0.0f, -6.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.5f, 0.5f, 0.5f});
	light.addComponent<DrawableComponent>(_renderer, shader, c.vertices, c.normals, c.indices, GL_TRIANGLES, RD_DEBUG);
	light.addComponent<LightComponent>(_renderer, glm::vec3{1.0f, 1.0f, 1.0f}, 0.3f);

	auto & light2(_ECS_manager->addEntity());
	light2.addComponent<TransformComponent>(glm::vec3{-7.0f, -2.0f, 6.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
	light2.addComponent<DrawableComponent>(_renderer, shader, c.vertices, c.normals, c.indices, GL_TRIANGLES, RD_DEBUG);
	light2.addComponent<LightComponent>(_renderer, glm::vec3{1.0f, 1.0f, 1.0f}, 0.2f);

	auto& mc(_ECS_manager->addEntity());
	mc.addComponent<TransformComponent>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
	mc.addComponent<DrawableComponent>(_renderer, shader, c.vertices, c.normals, c.indices, GL_TRIANGLES);
	mc.getComponent<DrawableComponent>().setColor({1.0f, 0.0f, 0.0f});
	auto& marchingCubeComponent = mc.addComponent<MarchingCubeComponent>(5.0f, 8.0f, 8.0f, 0.5f);

	auto& metaball(_ECS_manager->addEntity());
	metaball.addComponent<TransformComponent>(glm::vec3{0.0f, -1.0f, -1.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
	metaball.addComponent<MetaballComponent>(&marchingCubeComponent, 1.5f);
	auto& metaball2(_ECS_manager->addEntity());
	metaball2.addComponent<TransformComponent>(glm::vec3{0.0f, 1.0f, 1.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
	metaball2.addComponent<MetaballComponent>(&marchingCubeComponent, 1.5f);

	/*
	auto& metaball(_ECS_manager->addEntity());
	metaball.addComponent<TransformComponent>(glm::vec3{0.0f, -4.0f, -4.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
	metaball.addComponent<MetaballComponent>(&marchingCubeComponent, 1.5f);

	auto& metaball2(_ECS_manager->addEntity());
	metaball2.addComponent<TransformComponent>(glm::vec3{0.0f, -4.0f, 4.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
	metaball2.addComponent<MetaballComponent>(&marchingCubeComponent, 1.5f);

	auto& metaball3(_ECS_manager->addEntity());
	metaball3.addComponent<TransformComponent>(glm::vec3{0.0f, 4.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
	metaball3.addComponent<MetaballComponent>(&marchingCubeComponent, 1.7f);
	*/

}

void GlWidget::initializeGL()
{
	connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GlWidget::cleanup);
	if (!initializeOpenGLFunctions()) {
		exit(1);
	}
	glDebugMessageCallback(MessageCallback, 0);
	_renderer->initGl();
	_initScene();
}

void GlWidget::paintGL()
{
	_ECS_manager->update(_deltaTime);
	_renderer->draw(context()->defaultFramebufferObject());
	_endFrame();
	update();
}

void GlWidget::_endFrame()
{
	_printFPS();
	_computeDelta();
}

void GlWidget::_printFPS()
{
	std::uint64_t endTimerFPS = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	if (endTimerFPS-_startTimerFPS > 1000) {
		std::cout << "FPS: " << _frameCountSecond << std::endl;
		_frameCountSecond = 0;
		_startTimerFPS = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	}
	_frameCountSecond++;
}

void GlWidget::_computeDelta()
{
	if (_startTimerFrame != 0) {
		std::uint64_t endTimerFrame = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		_deltaTime = 1e-9*(endTimerFrame-_startTimerFrame);
	}
	_startTimerFrame = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
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
void GlWidget::wheelEvent(QWheelEvent *event) { _InputManager->wheelEvent(event); }

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

