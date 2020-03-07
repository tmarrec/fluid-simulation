#pragma once

#include "OpenGL.h"

class Light : public Entity {

public:
	Light(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
		std::string type, glm::vec3 color);
	~Light(void);
	
	void draw(glm::vec3 view_position, glm::mat4 projection, float delta_time,
				std::vector<std::shared_ptr<Entity>> lights) override;

	const std::string type() const;
	const glm::vec3 color() const;

private:
	glm::vec3 _color;
	const std::string _type;
};
