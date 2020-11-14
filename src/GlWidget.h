#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions_4_5_Core>
#include <QKeyEvent>
#include <chrono>
#include <cstdint>
#include <qevent.h>
#include <thread>

#include "utils.h"

#include "ECS.h"
#include "TransformComponent.h"
#include "CameraComponent.h"
#include "DrawableComponent.h"

class Renderer;
class ECS_Manager;
class InputManager;

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
	void keyReleaseEvent(QKeyEvent *event) override;
	void mousePressEvent(QMouseEvent *event) override;
	void mouseMoveEvent(QMouseEvent *event) override;
	void wheelEvent(QWheelEvent *event) override;

private:
	void _initScene();
	void _endFrame();
	void _computeDelta();
	void _printFPS();
	const std::shared_ptr<ECS_Manager> _ECS_manager;
	const std::shared_ptr<Renderer> _renderer;
	const std::shared_ptr<InputManager> _InputManager;

	std::uint64_t _frameCountSecond = 0;
	std::uint64_t _startTimerFPS = 0;
	std::uint64_t _startTimerFrame = 0;
	double _deltaTime = 0.0f;
};

