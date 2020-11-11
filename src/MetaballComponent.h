#pragma once

#include "ECS.h"
#include "TransformComponent.h"
#include "DrawableComponent.h"
#include "MarchingCube.h"
#include "glm/gtx/string_cast.hpp"
#include <cstdint>

class MetaballComponent : public Component
{
public:
	MetaballComponent(float __radius)
	: Component{}
	, _radius { __radius }
	{}

	void init() override
	{
		ASSERT(entity->hasComponent<TransformComponent>(), "entity should have a TransformComponent");
		ASSERT(entity->hasComponent<DrawableComponent>(), "entity should have a DrawableComponent");
		float cellSize = 1.0f;
		// Grid initialization
		for (float x = 0; x < 10; x += cellSize)
		{
			std::vector<std::vector<glm::vec3>> grid2D;
			for (float y = 0; y < 10; y += cellSize)
			{
				std::vector<glm::vec3> grid1D;
				for (float z = 0; z < 10; z += cellSize)
				{
					grid1D.emplace_back(glm::vec3{x,y,z});
				}
				grid2D.emplace_back(grid1D);
			}
			_grid.emplace_back(grid2D);
		}
		
		// TODO do this at update
		// check inside grid
		auto pos = entity->getComponent<TransformComponent>().position();
		for (std::uint64_t x = 0; x < _grid.size(); ++x)
		{
			for (std::uint64_t y = 0; y < _grid[x].size(); ++y)
			{
				for (std::uint64_t z = 0; z < _grid[x][y].size(); ++z)
				{
					float inside = std::pow(_radius, 2)/(std::pow(_grid[x][y][z].x-pos.x, 2)+std::pow(_grid[x][y][z].y-pos.y, 2)+std::pow(_grid[x][z][z].z-pos.z, 2));
					if (inside >= 1.0f)
					{
						std::cout << glm::to_string(_grid[x][y][z]) << std::endl;
					}
				}
			}
		}
	}
	void draw() override
	{
	}

	void update([[maybe_unused]] double _deltaTime) override
	{
		// TODO ici modifier vertices
		auto drawable = entity->getComponent<DrawableComponent>();
	}
	~MetaballComponent() override
	{
	}

private:
	float _radius;
	std::vector<std::vector<std::vector<glm::vec3>>> _grid;
};
