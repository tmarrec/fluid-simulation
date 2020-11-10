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
	DrawableComponent(Renderer__ __renderer, std::shared_ptr<Shader> __shader, std::vector<GLfloat> __vertices, std::vector<GLfloat> __normals, std::vector<GLuint> __indices, GLenum __drawMode, RDenum __debug=RD_NORMAL)
	: Component{}
	, _renderer { __renderer }
	, _vertices { std::make_shared<std::vector<GLfloat>>(__vertices) }
	, _normals { std::make_shared<std::vector<GLfloat>>(__normals) }
	, _indices { std::make_shared<std::vector<GLuint>>(__indices) }
	, _shader { __shader }
	, _drawMode { __drawMode }
	, _debug { __debug }
	{
		//TODO temp
		_color = {((double) rand() / (RAND_MAX)) + 1, ((double) rand() / (RAND_MAX)) + 1, ((double) rand() / (RAND_MAX)) + 1};
	}

	void init() override
	{
		_renderer->initDrawable(this);
	}
	void draw() override
	{
		_renderer->drawDrawable(this);	
	}
	void update([[maybe_unused]] double _deltaTime) override
	{}
	void updateGeometry()
	{
		_renderer->initDrawable(this);
	}

	void setVertices(std::vector<GLfloat> __vertices) {	_vertices = std::make_shared<std::vector<GLfloat>>(__vertices);	}
	void setNormals(std::vector<GLfloat> __normals)	{ _normals = std::make_shared<std::vector<GLfloat>>(__normals);	}
	void setIndices(std::vector<GLuint> __indices) { _indices = std::make_shared<std::vector<GLuint>>(__indices); }
	void setColor(glm::vec3 __color) { _color = __color; }

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
	RDenum debug() { return _debug; }
	glm::mat4 getModel()
	{
		ASSERT(entity->hasComponent<TransformComponent>(), "entity should have a TransformComponent");
		return entity->getComponent<TransformComponent>().getModel();
	}

	~DrawableComponent() override
	{
		_renderer->freeDrawable(this);
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

	glm::vec3 _color;
	std::shared_ptr<Shader> _shader;
	GLenum _drawMode;
	RDenum _debug;
};
