#include "src/CameraComponent.h"
#include "src/TransformComponent.h"
#include "utils.h"

#include "InputManager.h"
#include "ui/GlWidget.h"

InputManager::InputManager(GlWidget* __glWidget)
{
	_glWidget = __glWidget;	
	_input_timer = new QTimer(nullptr);
	connect(_input_timer, &QTimer::timeout, this, &InputManager::_process_inputs);
	_input_timer->start((1.0f/128.0f)*1000.0f);
}

void InputManager::_process_inputs()
{
	Entity* activeCameraEntity = nullptr;
	glm::vec3 cameraFrontVec;
	glm::vec3 cameraUpVec;
	glm::vec3 directionVec;
	if (_moveFront || _moveBack || _moveLeft || _moveRight || _moveUp || _moveDown)
	{
		activeCameraEntity = _glWidget->getActiveCamera();
		ASSERT(activeCameraEntity->hasComponent<TransformComponent>(), "activeCameraEntity should have a TransformComponent");
		ASSERT(activeCameraEntity->hasComponent<CameraComponent>(), "activeCameraEntity should have a CameraComponent");
		auto cameraComponent = activeCameraEntity->getComponent<CameraComponent>();
		cameraFrontVec = cameraComponent.front();
		cameraUpVec = cameraComponent.up();
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
		activeCameraEntity->getComponent<TransformComponent>().move(directionVec*activeCameraEntity->getComponent<CameraComponent>().speed());
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

void InputManager::mousePressEvent(QMouseEvent *event)
{
	_last_mouse_x = event->x();
	_last_mouse_y = event->y();
}

void InputManager::mouseMoveEvent(QMouseEvent *event)
{
	int x = event->x();
	int y = event->y();

	auto activeCameraEntity = _glWidget->getActiveCamera();
	ASSERT(activeCameraEntity->hasComponent<CameraComponent>(), "activeCameraEntity should have a CameraComponent");
	auto& activeCamera = activeCameraEntity->getComponent<CameraComponent>();
	float yaw = activeCamera.yaw();
	float pitch = activeCamera.pitch();

	float sensitivity = 0.3f;
	float offset_x = (x-_last_mouse_x)*sensitivity;
	float offset_y = (y-_last_mouse_y)*sensitivity;

	yaw += offset_x;
	pitch -= offset_y;
	pitch = glm::clamp(pitch, -89.9f, 89.9f);

	activeCamera.setYaw(yaw);
	activeCamera.setPitch(pitch);

	glm::vec3 dir;
	dir.x = cos(glm::radians(yaw))*cos(glm::radians(pitch));
	dir.y = sin(glm::radians(pitch));
	dir.z = sin(glm::radians(yaw))*cos(glm::radians(pitch));
	dir = glm::normalize(dir);
	activeCamera.setFront(dir);

	_last_mouse_x = x;
	_last_mouse_y = y;
}
