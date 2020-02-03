#include "GLWidget.h"
#include "shapes/Triangle.h"

#include <QMessageBox>
#include <QDateTime>

#include <iostream>

GLWidget::GLWidget(QWidget *parent)
		: QOpenGLWidget(parent)
{
	std::cout << "GLWidget constructor" << std::endl;

	setFixedSize(800, 600);
}

GLWidget::~GLWidget() {

}

void GLWidget::cleanup() {
	_openGL.reset(nullptr);
}

void GLWidget::mouseMoveEvent(QMouseEvent *event) {
	std::cout << "mouse move" << std::endl;
}

void GLWidget::initializeGL() {
	std::cout << "init_GL" << std::endl;
	connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GLWidget::cleanup);
	if (!initializeOpenGLFunctions()) {
		QMessageBox::critical(
						this,
						"OpenGL initialization error",
						"GLWidget::init_GL() : Unable to initialize OpenGL functions");
		exit(1);
	}	

	_openGL.reset(new Triangle{100, 150});
}

void GLWidget::paintGL() {
	std::int64_t start_time = QDateTime::currentMSecsSinceEpoch();
	_openGL->draw();
	glFinish();
	std::int64_t end_time = QDateTime::currentMSecsSinceEpoch();
	_last_time = end_time-start_time;
	std::cout << _last_time << std::endl;
}
