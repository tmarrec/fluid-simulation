#pragma once

class MainWindow;

#include <vector>
#include <memory>

#include "Entity.h"

class ECS {

public:
	ECS();
	~ECS(void);
	void add(std::unique_ptr<Entity> entity);
	void render_all(glm::vec3 view_position, glm::mat4 projection, float delta_time);
	void get(uint id);
	void move(uint id, char pos, uint value);

private:
	std::vector<std::unique_ptr<Entity>> _entities;
};
