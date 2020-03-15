#pragma once


class MainWindow;
class OpenGL;

#include "OpenGL.h"
#include "MainWindow.h"

#include <memory>

#include <QOpenGLWidget>
#include <QOpenGLFunctions_4_1_Core>


class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions_4_1_Core {
	Q_OBJECT

public:
	explicit GLWidget(MainWindow *parent = nullptr);	
	~GLWidget();

	void init();
	void cleanup();
	void make_current();
	void done_current();
	OpenGL * openGL();

protected:
	// Override from QOpenGLWidget
	void initializeGL() override;
	void paintGL() override;

private:
	OpenGL * _openGL;

	std::uint64_t _frame_count;
	std::uint64_t _start_timer_fps;
	std::uint64_t _start_timer_frame;
	
	MainWindow * _main_window;
};
