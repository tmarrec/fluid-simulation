#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Entity.h"

class Light : public Entity {

public:
	Light(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
		glm::vec3 color, MainWindow * main_window);
	~Light(void);
	
	void draw(glm::mat4 view, glm::mat4 projection, float delta_time,
				std::vector<std::shared_ptr<Entity>> lights) override;

	const glm::vec3 color() const;
	void set_shader(Shader shader) override;
	Shader & shader() override;
	Entity_Type type() override;

private:
	glm::vec3 _color;
};
