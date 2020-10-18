#pragma once

#include "ECS.h"

class TransformComponent : public Component
{
public:
	TransformComponent(glm::vec3 __position, glm::vec3 __rotation, glm::vec3 __scale)
	: Component{}
	, _position { __position }
	, _rotation { __rotation }
	, _scale { __scale }
	{}

	void move(glm::vec3 __moveVector)
	{
		_position += __moveVector;
	}
	void rotate(glm::vec3 __rotationVector)
	{
		_rotation += __rotationVector;
	}

	void setPosition(glm::vec3 __position) { _position = __position; }

	glm::vec3 position() const { return _position; }
	glm::vec3 rotation() const { return _rotation; }
	glm::vec3 scale() const { return _scale; }
	glm::mat4 getModel() const
	{
		glm::mat4 model {1.0f};
		model = glm::translate(model, _position);
		model = glm::rotate(model, glm::radians(_rotation.x), glm::vec3(1.0f, 0.0f, 0.0f));
		model = glm::rotate(model, glm::radians(_rotation.y), glm::vec3(0.0f, 1.0f, 0.0f));
		model = glm::rotate(model, glm::radians(_rotation.z), glm::vec3(0.0f, 0.0f, 1.0f));
		model = glm::scale(model, glm::vec3{_scale});
		return model;
	}

private:
	glm::vec3 _position;
	glm::vec3 _rotation;
	glm::vec3 _scale;
};
