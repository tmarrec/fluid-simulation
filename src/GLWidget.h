#pragma once


class MainWindow;
class OpenGL;

#include "OpenGL.h"
#include "MainWindow.h"

#include <memory>

#include <QOpenGLWidget>
#include <QOpenGLFunctions_4_5_Core>


class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions_4_5_Core {

public:
	explicit GLWidget(MainWindow *parent = nullptr);	
	~GLWidget();

	void cleanup();
	void add_shape(std::string shape);
	void select_entity(uint selected_id);
	void move(uint id, char pos, float value);

protected:
	// Override from QOpenGLWidget
	void initializeGL() override;
	void paintGL() override;

	void mouseMoveEvent(QMouseEvent *event) override;

private:
	std::unique_ptr<OpenGL> _openGL;

	std::uint64_t _frame_count;
	std::uint64_t _start_timer_fps;
	std::uint64_t _start_timer_frame;
	MainWindow * _main_window;
};
