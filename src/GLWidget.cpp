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
}

GLWidget::~GLWidget() {

}

void GLWidget::init() {
	setFixedSize(1280, 720);
	_openGL = new OpenGL{1280, 720, _main_window, this};
}

OpenGL * GLWidget::openGL() {
	return _openGL;
}

void GLWidget::cleanup() {

}

/*
void GLAPIENTRY
MessageCallback( GLenum source,
                 GLenum type,
                 GLuint id,
                 GLenum severity,
                 GLsizei length,
                 const GLchar* message,
                 const void* userParam )
{
	fprintf( stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
           ( type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : "" ),
            type, severity, message );
}
*/

void GLWidget::initializeGL() {
	connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GLWidget::cleanup);
	if (!initializeOpenGLFunctions()) {
		QMessageBox::critical(
						this,
						"OpenGL initialization error",
						"GLWidget::init_GL() : Unable to initialize OpenGL functions");
		exit(1);
	}	

	glEnable(GL_DEPTH_TEST);
	glViewport(0, 0, 1280, 720);
	//glEnable(GL_DEBUG_OUTPUT);
	//glDebugMessageCallback(MessageCallback, 0);

	std::cout << "QT version   : " << qVersion() << std::endl;
	std::cout << "Renderer     : " << glGetString(GL_RENDERER) << std::endl;
	std::cout << "Vendor       : " << glGetString(GL_VENDOR) << std::endl;
	std::cout << "Version      : " << glGetString(GL_VERSION) << std::endl;
	std::cout << "GLSL Version : " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
	std::cout << std::endl;

	_start_timer_fps = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

void GLWidget::make_current() {
	makeCurrent();
}

void GLWidget::done_current() {
	doneCurrent();
}

void GLWidget::paintGL() {
	_openGL->draw();
	glFinish();

	// Compte les FPS chaque secondes
	std::uint64_t end_timer_fps = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	if (end_timer_fps-_start_timer_fps > 1000) {
		_main_window->update_title_infos("FPS: " + std::to_string(_frame_count));
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
