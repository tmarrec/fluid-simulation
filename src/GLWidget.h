#pragma once

#include "OpenGL.h"

#include <memory>

#include <QOpenGLWidget>
#include <QOpenGLFunctions_4_5_Core>

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions_4_5_Core {

public:
	explicit GLWidget(QWidget *parent = nullptr);	
	~GLWidget();

	void cleanup();
	

protected:
	// Override from QOpenGLWidget
	void initializeGL() override;
	void paintGL() override;

	void mouseMoveEvent(QMouseEvent *event) override;

private:
	std::unique_ptr<OpenGL> _openGL;

	std::uint64_t _frame_count;
	std::uint64_t _start_timer;
};
