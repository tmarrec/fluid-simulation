#pragma once

#include "ECS.h"

class TransformComponent : public Component
{
private:
	glm::vec3 _position;
	glm::vec3 _rotation;
	glm::vec3 _scale;

public:
	TransformComponent(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
	: Component{}
	{
		_position = position;
		_rotation = rotation;
		_scale = scale;
	}

	glm::vec3 position()
	{
		return _position;
	}

	glm::vec3 rotation()
	{
		return _rotation;
	}

	glm::vec3 scale()
	{
		return _scale;
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

	~TransformComponent() override
	{
	}
};
