#include "GlWidget.h"
#include "src/ui/MainWindow.h"

#include "../Renderer.h"
#include "../TransformComponent.h"
#include "../CameraComponent.h"
#include "../DrawableComponent.h"

void GLAPIENTRY MessageCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam);

GlWidget::GlWidget(QWidget *parent)
: QOpenGLWidget{parent}
, _renderer{ std::shared_ptr<Renderer>(new Renderer) }
, _manager{ std::shared_ptr<Manager>(new Manager) }
{
	setFocus();
}

void GlWidget::_init()
{
	std::vector<GLfloat> v = std::vector<GLfloat>{
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
	};
	std::vector<GLfloat> n = std::vector<GLfloat>{
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
				};
	std::vector<GLuint> i = std::vector<GLuint>{
					0,  1,  2,  0,  2,  3,   //front
					4,  5,  6,  4,  6,  7,   //right
					8,  9,  10, 8,  10, 11,  //back
					12, 13, 14, 12, 14, 15,  //left
					16, 17, 18, 16, 18, 19,  //upper
					20, 21, 22, 20, 22, 23 	 //bottom
				};

	auto & cube(_manager->addEntity());
	cube.addComponent<TransformComponent>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{50.0f, 50.0f, 50.0f});
	cube.addComponent<DrawableComponent>(_renderer, "shaders/vert.vert", "shaders/frag.frag", std::make_shared<std::vector<GLfloat>>(v), std::make_shared<std::vector<GLfloat>>(n), std::make_shared<std::vector<GLuint>>(i));

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

	_init();
	
	//cout(std::string("QT Version     : ")+qVersion());
}

void GlWidget::paintGL()
{

	_renderer->setActiveCamera(&_camera->getComponent<CameraComponent>());
	_manager->update();
	_renderer->draw();


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
		case Qt::Key_R:
			{

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

