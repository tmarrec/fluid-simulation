#pragma once

#include "OpenGL.h"
#include "Shader.h"
#include "MainWindow.h"

class Shape : public Entity {

public:
	Shape(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
		std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> geometry,
		glm::vec3 color, Shader shader, MainWindow * main_window);
	~Shape(void);
	
	std::vector<GLuint> indices() const;

	void draw(glm::mat4 view, glm::mat4 projection, float delta_time,
			std::vector<std::shared_ptr<Entity>> lights) override;
	void draw_vertex();
	void use_shader(glm::mat4 view, glm::mat4 projection,
				std::vector<std::shared_ptr<Entity>> lights);

	void set_geometry();
	void set_shader(Shader shader) override;

	Entity_Type type() override;

	const glm::vec3 color() const;
	Shader & shader() override;
	
private:
	std::vector<GLfloat> 	_vertices;
	std::vector<GLfloat> 	_normals;
	std::vector<GLuint>		_indices;
	GLuint _VAO;

	glm::vec3 _color;
	Shader _shader;

};
