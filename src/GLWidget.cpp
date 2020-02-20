#include "GLWidget.h"

#include <QDateTime>
#include <QMessageBox>

#include <iostream>
#include <chrono>

GLWidget::GLWidget(MainWindow *parent)
			: QOpenGLWidget(parent)
			, _frame_count{0}
			, _start_timer_fps{0}
			, _start_timer_frame{0}
			, _main_window{parent}
{
	setFixedSize(1280, 720);
	_openGL = new OpenGL{1280, 720, _main_window};
}

GLWidget::~GLWidget() {

}

OpenGL * GLWidget::openGL() {
	return _openGL;
}

void GLWidget::cleanup() {

}

void GLWidget::move(uint id, char pos, float value) {
	_openGL->move(id, pos, value);
}

void GLWidget::initializeGL() {
	connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GLWidget::cleanup);
	if (!initializeOpenGLFunctions()) {
		QMessageBox::critical(
						this,
						"OpenGL initialization error",
						"GLWidget::init_GL() : Unable to initialize OpenGL functions");
		exit(1);
	}	

	_start_timer_fps = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

void GLWidget::paintGL() {
	_openGL->draw();
	glFinish();

	// Compte les FPS chaque secondes
	std::uint64_t end_timer_fps = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	if (end_timer_fps-_start_timer_fps > 1000) {
		std::cout << "FPS: " << _frame_count << std::endl;
		_frame_count = 0;
		_start_timer_fps = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	}

	// Calcule le Time Delta
	if (_start_timer_frame != 0) {
		std::uint64_t end_timer_frame = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		_openGL->set_delta_time(1e-9*(end_timer_frame-_start_timer_frame));
	}
	_start_timer_frame = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

	_frame_count++;
	update();
}
