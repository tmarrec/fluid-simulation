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


	std::int64_t _last_time;

};
