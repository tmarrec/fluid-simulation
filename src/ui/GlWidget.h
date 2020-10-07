#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions_4_5_Core>
#include <QKeyEvent>

#include "../utils.h"

class Renderer;
class Manager;
class Entity;

using Renderer__ = std::shared_ptr<Renderer>; 
using Manager__ = std::shared_ptr<Manager>; 

class GlWidget : public QOpenGLWidget, protected QOpenGLFunctions_4_5_Core
{
Q_OBJECT
public:
	explicit GlWidget(QWidget *parent = nullptr);

	~GlWidget();

public slots:
	void cleanup();

protected:
	void initializeGL() override;
	void paintGL() override;
	void resizeGL(int w, int h) override;

	void keyPressEvent(QKeyEvent *event) override;

private:
	void _init();
	const Renderer__ _renderer;
	const Manager__ _manager;
	Entity* _camera; //TODO Change that

};

