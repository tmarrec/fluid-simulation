#pragma once

#include <QKeyEvent>
#include <QTimer>
#include <QObject>
#include <qevent.h>
#include <memory>

class GlWidget;
class Renderer;

class InputManager : public QObject
{
Q_OBJECT
public:
	InputManager(GlWidget* __glWidget, std::shared_ptr<Renderer> __renderer);
	void keyPressEvent(QKeyEvent *event);
	void keyReleaseEvent(QKeyEvent *event);
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
	

private:
	void _process_inputs();

	bool _moveFront = false;
	bool _moveBack = false;
	bool _moveLeft = false;
	bool _moveRight = false;
	bool _moveUp = false;
	bool _moveDown = false;
	QTimer* _inputTimer;
	int _lastMouseX = 0;
	int _lastMouseY = 0;
	GlWidget* _glWidget;
	std::shared_ptr<Renderer> _renderer;
};
