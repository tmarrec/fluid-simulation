#pragma once

#include <GL/gl.h>

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
		std::shared_ptr<std::vector<GLfloat>> vertices, std::shared_ptr<std::vector<GLfloat>> normals,
		std::shared_ptr<std::vector<GLuint>> indices)
	: Component{}
	{
		_vertices = vertices;
		_normals = normals;
		_indices = indices;
		_renderer = __renderer;
		
		auto shader = new Shader{vertPath.c_str(), fragPath.c_str()};
		std::shared_ptr<Shader> uPtr {shader};
		_shader = uPtr;
	}

	void init() override
	{
		std::cout << "init" << std::endl;
		_renderer->initDrawable(this);
	}

	void draw() override
	{
		_renderer->useShader(this);	
	}

	void update() override
	{
	}

	~DrawableComponent() override
	{
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
