#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <string>
#include <vector>
#include <memory>

class Entity {

public:
	Entity(std::string name, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
	~Entity(void);

	virtual void draw(glm::vec3 view_position, glm::mat4 projection, float delta_time,
					std::vector<std::shared_ptr<Entity>> lights) = 0;

	void set_position(glm::vec3 position);
	void set_rotation(glm::vec3 rotation);
	glm::vec3 position() const;
	glm::vec3 rotation() const;

	void rotate_test(float delta_time);
	glm::mat4 get_model() const;
	glm::mat4 get_view(glm::vec3 view_position) const;

	unsigned long id() const;
	const std::string name() const;

private:
	glm::vec3 _position;
	glm::vec3 _rotation;
	glm::vec3 _scale;
	
	const unsigned long _id;
	static unsigned long _next_id;
	std::string _name;
};

