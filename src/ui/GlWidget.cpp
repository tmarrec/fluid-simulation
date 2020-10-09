#include "GlWidget.h"
#include "src/ui/MainWindow.h"

#include <cstdint>
#include <qevent.h>
#include <unistd.h>
#include <thread>

#include "../InputManager.h"
#include "../shapes.h"
#include "../Renderer.h"
#include "../TransformComponent.h"
#include "../CameraComponent.h"
#include "../DrawableComponent.h"

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

	Cube c;

	auto & cube(_ECS_manager->addEntity());
	cube.addComponent<TransformComponent>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	cube.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", c.vertices, c.normals, c.indices);

	auto & cube2(_ECS_manager->addEntity());
	cube2.addComponent<TransformComponent>(glm::vec3{20.0f, 10.0f, 80.0f}, glm::vec3{30.0f, 10.0f, 60.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	cube2.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", c.vertices, c.normals, c.indices);

	_activeCamera = &_ECS_manager->addEntity();
	_activeCamera->addComponent<CameraComponent>(0.0f, 0.0f, 5.0f, 90.0f);
	_activeCamera->addComponent<TransformComponent>(glm::vec3{-250.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f});

	_renderer->setActiveCamera(&_activeCamera->getComponent<CameraComponent>());
}

void GlWidget::initializeGL()
{
	connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GlWidget::cleanup);

	if (!initializeOpenGLFunctions()) {
		exit(1);
	}

	glEnable(GL_DEBUG_OUTPUT);
	glDebugMessageCallback(MessageCallback, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	
	_init();
	std::thread t ((Test(_ECS_manager)));
	t.detach();
}

void GlWidget::paintGL()
{
	_renderer->clear();
	_ECS_manager->update(_deltaTime);
	_ECS_manager->draw();
	glFinish();

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

void GlWidget::resizeGL(int w, int h) { _renderer->resizeGl(w, h); }
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

