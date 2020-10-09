#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions_4_5_Core>
#include <QKeyEvent>
#include <cstdint>

#include "../utils.h"

class Renderer;
class ECS_Manager;
class Entity;

using Renderer__ = std::shared_ptr<Renderer>; 
using ECS_Manager__ = std::shared_ptr<ECS_Manager>; 

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
	const ECS_Manager__ _manager;
	Entity* _camera; //TODO Change that

	std::uint64_t _frame_count = 0;
	std::uint64_t _start_timer_fps = 0;
	std::uint64_t _start_timer_frame = 0;

};

