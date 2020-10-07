#include "GlWidget.h"
#include "src/ui/MainWindow.h"
#include <memory>


GlWidget::GlWidget(QWidget *parent)
: QOpenGLWidget{parent}
{

}

void GLAPIENTRY
MessageCallback( GLenum source,
                 GLenum type,
                 GLuint id,
                 GLenum severity,
                 GLsizei length,
                 const GLchar* message,
                 const void* userParam )
{
	if (severity != 33387)
	{
		fprintf(stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
			(type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""), type, severity, message);
	}
}

void GlWidget::initializeGL()
{
	connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GlWidget::cleanup);

	if (!initializeOpenGLFunctions()) {
		exit(1);
	}

	glEnable              ( GL_DEBUG_OUTPUT );
	glDebugMessageCallback( MessageCallback, 0 );

	//cout(std::string("QT Version     : ")+qVersion());
}

void GlWidget::paintGL()
{
	std::cout << "paint" << std::endl;
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

void GlWidget::setRenderer(Renderer__ __renderer)
{
	ASSERT(__renderer, "__renderer should not be nullptr");
	_renderer = __renderer;
}

GlWidget::~GlWidget()
{

}

