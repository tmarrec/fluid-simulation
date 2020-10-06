#pragma once

#include "ECS.h"

#include "TransformComponent.h"
#include <GL/gl.h>

class DrawableComponent : public Component
{
private:
	std::shared_ptr<std::vector<GLfloat>> _vertices;
	std::shared_ptr<std::vector<GLfloat>> _normals;
	std::shared_ptr<std::vector<GLuint>>  _indices;

	GLuint _VAO;
	GLuint _VBO;
	GLuint _NBO;
	GLuint _EBO;

	glm::vec3 _color = {1.0f, 0.5f, 0.5f};
	std::shared_ptr<Shader> _shader;


public:
	DrawableComponent(std::string vertPath, std::string fragPath,
		std::shared_ptr<std::vector<GLfloat>> vertices, std::shared_ptr<std::vector<GLfloat>> normals,
		std::shared_ptr<std::vector<GLuint>> indices)
	: Component{}
	{
		_vertices = vertices;
		_normals = normals;
		_indices = indices;
		
		auto shader = new Shader{vertPath.c_str(), fragPath.c_str()};
		std::shared_ptr<Shader> uPtr {shader};
		_shader = uPtr;
	}

	void init() override
	{
		// TODO maybe should send pointers
		Message initMsg {INIT_DRAWABLE, _vertices, _normals, _indices, _VAO,
			_VBO, _NBO, _EBO, _color, _shader};
		entity->entityPostMessage(initMsg);
	}

	void draw() override
	{
		if (!entity->hasComponent<TransformComponent>())
		{
			std::cout << "ERROR: DrawableComponent.h : The entity does not have a TransformComponent" << std::endl;
			exit(1);
		}

		auto position = entity->getComponent<TransformComponent>().position();
		auto rotation = entity->getComponent<TransformComponent>().rotation();
		auto scale = entity->getComponent<TransformComponent>().scale();

		// TODO maybe should send pointers
		Message drawMsg {ASK_CAMERA_INFOS_FOR_DRAW, _vertices, _normals, _indices, _VAO,
			_VBO, _NBO, _EBO, _color, _shader, position, rotation, scale};
		entity->entityPostMessage(drawMsg);
	}

	void update() override
	{
	}

	~DrawableComponent() override
	{
		Message freeMsg {FREE_DRAWABLE, _vertices, _normals, _indices, _VAO,
			_VBO, _NBO, _EBO, _color, _shader};
		entity->entityPostMessage(freeMsg);
	}

	std::shared_ptr<std::vector<GLfloat>> vertices()
	{
		return _vertices;
	}

	std::shared_ptr<std::vector<GLfloat>> normals()
	{
		return _normals;
	}

	std::shared_ptr<std::vector<GLuint>> indices()
	{
		return _indices;
	}

	glm::vec3 color() const
	{
		return _color;
	}

	std::shared_ptr<Shader> shader() const
	{
		return _shader;
	}
};
