#include "utils.h"

#include "InputManager.h"
#include "ui/GlWidget.h"
#include <qwindowdefs.h>

InputManager::InputManager(GlWidget* __glWidget)
{
	_glWidget = __glWidget;	
	_input_timer = new QTimer(nullptr);
	connect(_input_timer, &QTimer::timeout, this, &InputManager::_process_inputs);
	_input_timer->start((1.0f/128.0f)*1000.0f);

}

void InputManager::_process_inputs()
{
	if (_moveFront)
	{
		glm::vec3 directionVec = {1.0f, 0.0f, 0.0f};
		auto activeCamera = _glWidget->getActiveCamera();
		activeCamera->getComponent<TransformComponent>().move(directionVec*activeCamera->getComponent<CameraComponent>().speed());
	}
	if (_moveBack)
	{
		glm::vec3 directionVec = {-1.0f, 0.0f, 0.0f};
		auto activeCamera = _glWidget->getActiveCamera();
		activeCamera->getComponent<TransformComponent>().move(directionVec*activeCamera->getComponent<CameraComponent>().speed());
	}
	if (_moveLeft)
	{
		glm::vec3 directionVec = {0.0f, 0.0f, -1.0f};
		auto activeCamera = _glWidget->getActiveCamera();
		activeCamera->getComponent<TransformComponent>().move(directionVec*activeCamera->getComponent<CameraComponent>().speed());
	}
	if (_moveRight)
	{
		glm::vec3 directionVec = {0.0f, 0.0f, 1.0f};
		auto activeCamera = _glWidget->getActiveCamera();
		activeCamera->getComponent<TransformComponent>().move(directionVec*activeCamera->getComponent<CameraComponent>().speed());
	}
	if (_moveUp)
	{
		glm::vec3 directionVec = {0.0f, 1.0f, 0.0f};
		auto activeCamera = _glWidget->getActiveCamera();
		activeCamera->getComponent<TransformComponent>().move(directionVec*activeCamera->getComponent<CameraComponent>().speed());
	}
	if (_moveDown)
	{
		glm::vec3 directionVec = {0.0f, -1.0f, 0.0f};
		auto activeCamera = _glWidget->getActiveCamera();
		activeCamera->getComponent<TransformComponent>().move(directionVec*activeCamera->getComponent<CameraComponent>().speed());
	}
}

void InputManager::keyPressEvent(QKeyEvent *event)
{
	switch(event->key())
	{
		case Qt::Key_Z:
			_moveFront = true;
			break;
		case Qt::Key_S:
			_moveBack = true;
			break;
		case Qt::Key_Q:
			_moveLeft = true;
			break;
		case Qt::Key_D:
			_moveRight = true;
			break;
		case Qt::Key_Shift:
			_moveUp = true;
			break;
		case Qt::Key_Control:
			_moveDown = true;
			break;

		default:
			break;
	}
}

void InputManager::keyReleaseEvent(QKeyEvent *event)
{
	switch(event->key())
	{
		case Qt::Key_Z:
			_moveFront = false;
			break;
		case Qt::Key_S:
			_moveBack = false;
			break;
		case Qt::Key_Q:
			_moveLeft = false;
			break;
		case Qt::Key_D:
			_moveRight = false;
			break;
		case Qt::Key_Shift:
			_moveUp = false;
			break;
		case Qt::Key_Control:
			_moveDown = false;
			break;

		default:
			break;
	}
}

