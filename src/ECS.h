#pragma once

#include <vector>
#include <memory>

#include "Entity.h"

class ECS {

public:
	ECS(void);
	~ECS(void);
	void add(std::unique_ptr<Entity> entity);
	void render_all(glm::vec3 view_position, glm::mat4 projection, float delta_time);

private:
	std::vector<std::unique_ptr<Entity>> _entities;
};
