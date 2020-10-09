#include "GlWidget.h"
#include "src/ui/MainWindow.h"

#include "../shapes.h"
#include "../Renderer.h"
#include "../TransformComponent.h"
#include "../CameraComponent.h"
#include "../DrawableComponent.h"

void GLAPIENTRY MessageCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam);

GlWidget::GlWidget(QWidget *parent)
: QOpenGLWidget{parent}
, _renderer{ std::shared_ptr<Renderer>(new Renderer) }
, _manager{ std::shared_ptr<ECS_Manager>(new ECS_Manager) }
{
	setFocus();
}

void GlWidget::_init()
{

	Cube c;

	auto & cube(_manager->addEntity());
	cube.addComponent<TransformComponent>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	cube.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", c.vertices, c.normals, c.indices);

	auto & cube2(_manager->addEntity());
	cube2.addComponent<TransformComponent>(glm::vec3{20.0f, 10.0f, 80.0f}, glm::vec3{30.0f, 10.0f, 60.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	cube2.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", c.vertices, c.normals, c.indices);

	_camera = &_manager->addEntity();
	_camera->addComponent<CameraComponent>(0.0f, 0.0f, 15.0f, 90.0f);
	_camera->addComponent<TransformComponent>(glm::vec3{-250.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
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
	
	//cout(std::string("QT Version     : ")+qVersion());
}

void GlWidget::paintGL()
{

	_manager->update(); //TODO a faire dans la gameloop

	_renderer->setActiveCamera(&_camera->getComponent<CameraComponent>());
	_renderer->clear();
	_manager->draw();

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
		//_openGL->set_delta_time(1e-9*(end_timer_frame-_start_timer_frame));
	}
	_start_timer_frame = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	_frame_count++;
	//#####################################

	update();
}

void GlWidget::resizeGL(int w, int h)
{
	_renderer->resizeGl(w, h);
}

void GlWidget::keyPressEvent(QKeyEvent *event)
{
	switch(event->key())
	{
		case Qt::Key_Z:
			{
				glm::vec3 directionVec = {1.0f, 0.0f, 0.0f};
				_camera->getComponent<TransformComponent>().move(directionVec*_camera->getComponent<CameraComponent>().speed());
			}
			break;

		default:
			break;
	}
}

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

