#include "GlWidget.h"


GlWidget::GlWidget(QWidget *parent)
: QOpenGLWidget{parent}
, System{qobject_cast<MainWindow*>(parentWidget()->parentWidget())->__messageBus()}
, _glWidth{0}
, _glHeight{0}
{
	Message helloMsg {HELLO, this};
	postMessage(helloMsg);
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
	fprintf(stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
		(type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""), type, severity, message);
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
	//Message initGLMsg {INIT_GL, size().width(), size().height(), context()};
	Message initGLMsg {INIT_GL, size().width(), size().height(), this};
	postMessage(initGLMsg);

}

void GlWidget::paintGL()
{
	Message drawMsg {ASK_ENTITIES_DRAW};
	postMessage(drawMsg);
}

void GlWidget::resizeGL(int w, int h)
{
	_glWidth = w;
	_glHeight = h;
	Message resizeMsg {RESIZE_GL, w, h};
	postMessage(resizeMsg);
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

void GlWidget::handleMessage(Message & msg)
{
	switch(msg._type)
	{
		case HELLO_ACK:
			cout("Loaded in the \033[45m\033[1m[MessageBus]\033[49m\033[0m");
			break;

		case ASK_GLWIDGET_SIZE:
			{
				Message glSizeMsg {GL_SIZE, _glWidth, _glHeight};
				postMessage(glSizeMsg);
			}
			break;

		default:
			break;
	}
}

void GlWidget::cout(std::string string) const
{
	std::cout << "0x" << std::hex << std::this_thread::get_id() << " ";
	std::cout << "  \033[46m\033[1m";
	std::cout << "[GlWidget]";
	std::cout << "\033[49m\033[0m";
	std::cout << " " << string << std::endl;

	std::cout << context() << std::endl;
}

void GlWidget::cleanup()
{

}

GlWidget::~GlWidget()
{

}

