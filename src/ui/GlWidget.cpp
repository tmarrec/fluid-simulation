#include "GlWidget.h"

#include <iostream>
#include "../MainWindow.h"

GlWidget::GlWidget(QWidget *parent)
: QOpenGLWidget{parent}
, System{qobject_cast<MainWindow*>(parentWidget()->parentWidget())->__messageBus()}
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
}

void GlWidget::paintGL()
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glFinish();
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

