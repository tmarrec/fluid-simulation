#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions_4_5_Core>
#include <QKeyEvent>

#include "../Renderer.h"

using Renderer__ = std::shared_ptr<Renderer>; 

class GlWidget : public QOpenGLWidget, protected QOpenGLFunctions_4_5_Core
{
Q_OBJECT
public:
	explicit GlWidget(QWidget *parent = nullptr);

	void setRenderer(Renderer__ __renderer);

	~GlWidget();

public slots:
	void cleanup();

protected:
	void initializeGL() override;
	void paintGL() override;
	void resizeGL(int w, int h) override;

	void keyPressEvent(QKeyEvent *event) override;

private:
	Renderer__ _renderer = nullptr;

};

