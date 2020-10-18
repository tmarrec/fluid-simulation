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
public:
	DrawableComponent(Renderer__ __renderer, std::string __vertPath, std::string __fragPath, std::vector<GLfloat> __vertices, std::vector<GLfloat> __normals, std::vector<GLuint> __indices, GLenum __drawMode)
	: Component{}
	, _renderer { __renderer }
	, _vertices { std::make_shared<std::vector<GLfloat>>(__vertices) }
	, _normals { std::make_shared<std::vector<GLfloat>>(__normals) }
	, _indices { std::make_shared<std::vector<GLuint>>(__indices) }
	, _shader { std::make_shared<Shader>(Shader{__vertPath.c_str(), __fragPath.c_str()}) }
	, _drawMode { __drawMode }
	{}

	void init() override
	{
		_renderer->initDrawable(this);
	}
	void draw() override
	{
		_renderer->drawDrawable(this);	
	}
	void update([[maybe_unused]] double _deltaTime) override
	{
		//entity->getComponent<TransformComponent>().rotate({0.0f, 0.01f, 0.0f});
		//auto test = glm::vec3{0.0f, 0.0f, 100.0f}*(float)_deltaTime;
		//entity->getComponent<TransformComponent>().move(test);
	}
	~DrawableComponent() override
	{
		_renderer->freeDrawable(this);
	}

	std::shared_ptr<std::vector<GLfloat>> vertices() const { return _vertices; }
	std::shared_ptr<std::vector<GLfloat>> normals() const { return _normals; }
	std::shared_ptr<std::vector<GLuint>> indices() const { return _indices; }
	glm::vec3 color() const	{ return _color; }
	std::shared_ptr<Shader> shader() const { return _shader; }
	GLenum drawMode() const { return _drawMode; }
	GLuint& VAO() { return _VAO; }
	GLuint& VBO() { return _VBO; }
	GLuint& NBO() { return _NBO; }
	GLuint& EBO() { return _EBO; }
	glm::mat4 getModel()
	{
		ASSERT(entity->hasComponent<TransformComponent>(), "entity should have a TransformComponent");
		return entity->getComponent<TransformComponent>().getModel();
	}

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
	GLenum _drawMode;
};
