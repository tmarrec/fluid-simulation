#include "GLWidget.h"
#include "shapes/Triangle.h"
#include "MainWindow.h"

#include <QDateTime>
#include <QMessageBox>

#include <iostream>
#include <chrono>

GLWidget::GLWidget(QWidget *parent)
			: QOpenGLWidget(parent)
			, _frame_count{0}
			, _start_timer{0}
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

	_start_timer = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

	_openGL.reset(new Triangle{100, 150});
}

void GLWidget::paintGL() {
	_openGL->draw();
	glFinish();
	std::uint64_t end_timer = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

	if (end_timer-_start_timer > 1000) {
		std::cout << "FPS: " << _frame_count << std::endl;
		_frame_count = 0;
		_start_timer = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	}

	_frame_count++;
	update();
}