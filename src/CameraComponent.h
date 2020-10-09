#pragma once

#include "ECS.h"

#include "TransformComponent.h"

class CameraComponent : public Component
{
private:
	float _yaw;
	float _pitch;
	float _speed;
	float _FOV;
	glm::vec3 _front = {1.0f, 0.0f, 0.0f};
	glm::vec3 _up = {0.0f, 1.0f, 0.0f};
	glm::mat4 _projection;

public:
	CameraComponent(float yaw, float pitch, float speed,
		float FOV)
	: Component{}
	{
		_yaw = yaw;
		_pitch = pitch;
		_speed = speed;
		_FOV = FOV;
		_projection = glm::mat4{1.0f};
	}

	float yaw() const
	{
		return _yaw;
	}
	
	float pitch() const
	{
		return _pitch;
	}

	float speed() const
	{
		return _speed;
	}
	
	glm::vec3 front() const
	{
		return _front;
	}

	glm::vec3 up() const
	{
		return _up;
	}

	glm::mat4 view() const
	{
		ASSERT(entity->hasComponent<TransformComponent>(), "entity should have a TransformComponent");
		auto position = entity->getComponent<TransformComponent>().position();
		return glm::lookAt(position, position+_front, _up);	
	}

	glm::mat4 projection() const
	{
		return _projection;
	}

	void setProjection(int width, int height)
	{
		_projection = glm::infinitePerspective(glm::radians(_FOV), (float)width/(float)height, 0.1f);
	}

	void setYaw(float __yaw)
	{
		_yaw = __yaw;
	}

	void setPitch(float __pitch)
	{
		_pitch = __pitch;
	}

	void setFront(glm::vec3 __front)
	{
		_front = __front;
	}

	void init() override
	{
	}

	void draw() override
	{
	}

	void update(double _deltaTime) override
	{
	}

	~CameraComponent() override
	{
	}
};
