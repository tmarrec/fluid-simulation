#pragma once

class MainWindow;

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Shader.h"
#include <string>
#include <vector>
#include <memory>

enum Entity_Type { SHAPE, CAMERA, LIGHT };


class Entity {

public:
	Entity(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window);
	~Entity(void);

	virtual void draw(glm::mat4 view_position, glm::mat4 projection, float delta_time,
					std::vector<std::shared_ptr<Entity>> lights) = 0;

	virtual void set_shader(Shader shader) = 0; 
	virtual Shader & shader() = 0;
	virtual Entity_Type type() = 0;

	void set_position(glm::vec3 position);
	void set_rotation(glm::vec3 rotation);
	void set_scale(glm::vec3 scale);
	glm::vec3 position() const;
	glm::vec3 rotation() const;
	glm::vec3 scale() const;

	void rotate_test(float delta_time);
	glm::mat4 get_model() const;
	glm::mat4 get_view(glm::vec3 view_position) const;

	unsigned long id() const;

private:
	glm::vec3 _position;
	glm::vec3 _rotation;
	glm::vec3 _scale;
	
	const unsigned long _id;
	static unsigned long _next_id;
	MainWindow * _main_window;
};

