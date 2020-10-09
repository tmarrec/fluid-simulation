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

	void move(glm::vec3 moveVector)
	{
		_position += moveVector;
	}

	glm::mat4 getModel()
	{
		glm::mat4 model {1.0f};
		model = glm::translate(model, _position);
		model = glm::rotate(model, glm::radians(_rotation.x), glm::vec3(1.0f, 0.0f, 0.0f));
		model = glm::rotate(model, glm::radians(_rotation.y), glm::vec3(0.0f, 1.0f, 0.0f));
		model = glm::rotate(model, glm::radians(_rotation.z), glm::vec3(0.0f, 0.0f, 1.0f));
		model = glm::scale(model, glm::vec3{_scale});

		return model;
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

	~TransformComponent() override
	{
	}
};
