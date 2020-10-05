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
		// need width and height :(
		//_projection = glm::perspective(glm::radians(FOV), (float)
	}

	glm::mat4 view() const
	{
		if (entity->hasComponent<TransformComponent>())
		{
			auto position = entity->getComponent<TransformComponent>().position();
			return glm::lookAt(position, position+_front, _up);	
		}
		else
		{
			std::cout << "ERROR: CameraComponent.h : The entity does not have a TransformComponent" << std::endl;
			exit(1);
		}
	}

	void init() override
	{
	}

	void draw() override
	{
	}

	void update() override
	{
	}

	~CameraComponent() override
	{
	}
};
