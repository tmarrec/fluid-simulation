#pragma once

#include <vector>

#include "shapes/Triangle.h"

class ECS {

public:
	ECS(void);
	~ECS(void);

private:
	std::vector<std::unique_ptr<Triangle>> _entities;
};
