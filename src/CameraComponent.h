#pragma once

#include "ECS.h"

#include "TransformComponent.h"

class CameraComponent : public Component
{
public:
	CameraComponent(float __yaw, float __pitch, float __speed, float __FOV)
	: Component{}
	, _yaw { __yaw }
	, _pitch { __pitch }
	, _speed { __speed }
	, _FOV { __FOV }
	, _projection { glm::mat4{1.0f} }
	{}

	float yaw() const { return _yaw; }
	float pitch() const { return _pitch;	}
	float speed() const	{ return _speed; }
	glm::vec3 front() const { return _front; }
	glm::vec3 up() const { return _up; }
	glm::mat4 projection() const { return _projection; }
	glm::mat4 getView() const {
		ASSERT(entity->hasComponent<TransformComponent>(), "entity should have a TransformComponent");
		auto position = entity->getComponent<TransformComponent>().position();
		return glm::lookAt(position, position+_front, _up);	
	}

	void setYaw(float __yaw) { _yaw = __yaw; }
	void setPitch(float __pitch) { _pitch = __pitch; }
	void setFront(glm::vec3 __front) { _front = __front; }
	void setProjection(int __width, int __height)
	{
		_projection = glm::infinitePerspective(glm::radians(_FOV), (float)__width/(float)__height, 0.1f);
	}

private:
	float _yaw;
	float _pitch;
	float _speed;
	float _FOV;
	glm::vec3 _front = {1.0f, 0.0f, 0.0f};
	glm::vec3 _up = {0.0f, 1.0f, 0.0f};
	glm::mat4 _projection;
};
