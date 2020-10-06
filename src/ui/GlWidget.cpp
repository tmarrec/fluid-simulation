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

void GlWidget::initializeGL()
{
	connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GlWidget::cleanup);

	if (!initializeOpenGLFunctions()) {
		exit(1);
	}

	//cout(std::string("QT Version     : ")+qVersion());
	Message initGLMsg {INIT_GL, size().width(), size().height()};
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
				Message test {TEST};
				postMessage(test);
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
}

void GlWidget::cleanup()
{

}

GlWidget::~GlWidget()
{

}

