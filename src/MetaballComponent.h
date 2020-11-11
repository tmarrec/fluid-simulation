#pragma once

#include "ECS.h"
#include "TransformComponent.h"
#include "MarchingCube.h"
#include <cstdint>

class MetaballComponent : public Component
{
public:
	MetaballComponent(MarchingCubeComponent* __marchingCubeComponent, float __radius)
	: Component{}
	, _marchingCubeComponent { __marchingCubeComponent }
	, _radius { __radius }
	{}

	void init() override
	{
		ASSERT(entity->hasComponent<TransformComponent>(), "entity should have a TransformComponent");
	}

	void update([[maybe_unused]] double _deltaTime) override
	{
		if (entity->getEntityID() == 3)
		{
			auto& transformComponent = entity->getComponent<TransformComponent>();
			//transformComponent.move({0.0f, 0.0f, 0.01f});
			//_needUpdate = true;
		}
		
		auto pos = entity->getComponent<TransformComponent>().position();
		if (!_needUpdate) return;
		_needUpdate = false;

		float isoLevel = 1.0f;
		auto& grid = _marchingCubeComponent->grid();
		for (std::uint64_t x = 0; x < grid.size(); ++x)
		{
			for (std::uint64_t y = 0; y < grid[x].size(); ++y)
			{
				for (std::uint64_t z = 0; z < grid[x][y].size(); ++z)
				{
					// loop cell points	
					for (std::uint64_t i = 0; i < 8; ++i)
					{
						float inside = std::pow(_radius, 2)/(std::pow(grid[x][y][z].points[i].x-pos.x, 2)+std::pow(grid[x][y][z].points[i].y-pos.y, 2)+std::pow(grid[x][y][z].points[i].z-pos.z, 2));
						if (inside >= 1.0f)
						{
							glm::vec3 center = pos;
							_marchingCubeComponent->changeGrid(x, y, z, i, isoLevel, center);
						}
					}
				}
			}
		}
	}

	~MetaballComponent() override
	{
	}

private:
	MarchingCubeComponent* _marchingCubeComponent;
	float _radius;
	bool _needUpdate = true;
};
