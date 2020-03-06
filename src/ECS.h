#pragma once

#include <vector>
#include <memory>

#include "Entity.h"

/*
class ECS : public QObject {
	Q_OBJECT
*/
class ECS {

public:
	ECS();
	~ECS(void);
	void add(std::shared_ptr<Entity> entity);
	void render_all(glm::vec3 view_position, glm::mat4 projection, float delta_time);

private:
	std::vector<std::shared_ptr<Entity>> _entities;
};
