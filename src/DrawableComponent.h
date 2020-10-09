#pragma once

#include <GL/gl.h>
#include <memory>

#include "ECS.h"
#include "TransformComponent.h"
#include "Shader.h"
#include "Renderer.h"

using Renderer__ = std::shared_ptr<Renderer>; 

class DrawableComponent : public Component
{
private:
	Renderer__ _renderer;
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
	DrawableComponent(Renderer__ __renderer, std::string vertPath, std::string fragPath,
		std::vector<GLfloat> vertices, std::vector<GLfloat> normals,
		std::vector<GLuint> indices)
	: Component{}
	{
		_vertices = std::make_shared<std::vector<GLfloat>>(vertices);
		_normals = std::make_shared<std::vector<GLfloat>>(normals);
		_indices = std::make_shared<std::vector<GLuint>>(indices);
		_renderer = __renderer;
		
		auto shader = new Shader{vertPath.c_str(), fragPath.c_str()};
		std::shared_ptr<Shader> uPtr {shader};
		_shader = uPtr;
	}

	void init() override
	{
		_renderer->initDrawable(this);
	}

	void draw() override
	{
		_renderer->drawDrawable(this);	
	}

	void update() override
	{
	}

	~DrawableComponent() override
	{
		_renderer->freeDrawable(this);
	}

	glm::mat4 getModel()
	{
		ASSERT(entity->hasComponent<TransformComponent>(), "entity should have a TransformComponent");
		return entity->getComponent<TransformComponent>().getModel();
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

	GLuint& VAO()
	{
		return _VAO;
	}

	GLuint& VBO()
	{
		return _VBO;
	}

	GLuint& NBO()
	{
		return _NBO;
	}

	GLuint& EBO()
	{
		return _EBO;
	}
};
