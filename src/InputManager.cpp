#include "src/CameraComponent.h"
#include "src/TransformComponent.h"
#include "utils.h"

#include "InputManager.h"
#include "GlWidget.h"
#include "Renderer.h"
#include <qnamespace.h>

InputManager::InputManager(GlWidget* __glWidget, std::shared_ptr<Renderer> __renderer)
: _glWidget { __glWidget }
, _renderer { __renderer }
{
	_inputTimer = new QTimer(nullptr);
	connect(_inputTimer, &QTimer::timeout, this, &InputManager::_process_inputs);
	_inputTimer->start((1.0f/128.0f)*1000.0f);
}

void InputManager::_process_inputs()
{
	CameraComponent* cameraComponent = _renderer->getActiveCamera();
	ASSERT(cameraComponent, "cameraComponent should not be nullptr");
	glm::vec3 cameraFrontVec;
	glm::vec3 cameraUpVec;
	glm::vec3 directionVec;
	if (_moveFront || _moveBack || _moveLeft || _moveRight || _moveUp || _moveDown)
	{
		cameraFrontVec = cameraComponent->front();
		cameraUpVec = cameraComponent->up();
	}

	if (_moveFront)
	{
		directionVec = glm::normalize(cameraFrontVec);
	}
	if (_moveBack)
	{
		directionVec = -glm::normalize(cameraFrontVec);
	}
	if (_moveLeft)
	{
		directionVec = -glm::normalize(glm::cross(cameraFrontVec, cameraUpVec));
	}
	if (_moveRight)
	{
		directionVec = glm::normalize(glm::cross(cameraFrontVec, cameraUpVec));
	}
	if (_moveUp)
	{
		directionVec = {0.0f, 1.0f, 0.0f};
	}
	if (_moveDown)
	{
		directionVec = {0.0f, -1.0f, 0.0f};
	}

	if (_moveFront || _moveBack || _moveLeft || _moveRight || _moveUp || _moveDown)
	{
		cameraComponent->entity->getComponent<TransformComponent>().move(directionVec*cameraComponent->entity->getComponent<CameraComponent>().speed());
	}
}

void InputManager::wheelEvent(QWheelEvent *event)
{
	CameraComponent* cameraComponent = _renderer->getActiveCamera();
	ASSERT(cameraComponent, "cameraComponent should not be nullptr");
	if (event->angleDelta().y() > 0)
	{
		cameraComponent->changeFOV(-5.0f);
	}
	else
	{
		cameraComponent->changeFOV(5.0f);
	}
}

void InputManager::keyPressEvent(QKeyEvent *event)
{
	_glWidget->makeCurrent();
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
		case Qt::Key_P:
			_renderer->switchPolygonmode();
			break;
		case Qt::Key_N:
			_renderer->switchShowNormals();
		default:
			break;
	}
	_glWidget->doneCurrent();
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

void InputManager::mousePressEvent(QMouseEvent *event)
{
	_lastMouseX = event->x();
	_lastMouseY = event->y();
}

void InputManager::mouseMoveEvent(QMouseEvent *event)
{
	int x = event->x();
	int y = event->y();

	CameraComponent* cameraComponent = _renderer->getActiveCamera();
	ASSERT(cameraComponent, "cameraComponent should not be nullptr");
	float yaw = cameraComponent->yaw();
	float pitch = cameraComponent->pitch();

	float sensitivity = 0.3f;
	float offset_x = (x-_lastMouseX)*sensitivity;
	float offset_y = (y-_lastMouseY)*sensitivity;

	yaw += offset_x;
	pitch -= offset_y;
	pitch = glm::clamp(pitch, -89.9f, 89.9f);

	cameraComponent->setYaw(yaw);
	cameraComponent->setPitch(pitch);

	glm::vec3 dir;
	dir.x = cos(glm::radians(yaw))*cos(glm::radians(pitch));
	dir.y = sin(glm::radians(pitch));
	dir.z = sin(glm::radians(yaw))*cos(glm::radians(pitch));
	dir = glm::normalize(dir);
	cameraComponent->setFront(dir);

	_lastMouseX = x;
	_lastMouseY = y;
}
