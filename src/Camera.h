#pragma once

#include <string>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Entity.h"

class Camera : public Entity {

public:
	Camera(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
		std::string type, float FOV);
	~Camera(void);
	
	void draw(glm::vec3 view_position, glm::mat4 projection, float delta_time,
				std::vector<std::shared_ptr<Entity>> lights) override;
	
	const std::string type() const;
	float FOV() const;

private:
	const std::string _type;
	float _FOV;
};
