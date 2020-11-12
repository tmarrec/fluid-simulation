#pragma once

#include "ECS.h"
#include "TransformComponent.h"
#include "MarchingCubeComponent.h"
#include <cstdint>
#include <chrono>

class MetaballComponent : public Component
{
public:
	MetaballComponent(MarchingCubeComponent* __marchingCubeComponent, float __radius)
	: Component{}
	, _marchingCubeComponent { __marchingCubeComponent }
	, _radius { __radius }
	{}

	void update([[maybe_unused]] double _deltaTime) override
	{
		if (entity->getEntityID() == 3)
		{
			auto& transformComponent = entity->getComponent<TransformComponent>();
			transformComponent.move({0.0f, 0.0f, 0.005f});
		}
		else 
		{
			auto& transformComponent = entity->getComponent<TransformComponent>();
			transformComponent.move({0.0f, 0.0f, -0.005f});
		}
		std::uint64_t start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		ASSERT(entity->hasComponent<TransformComponent>(), "entity should have a TransformComponent");
		auto pos = entity->getComponent<TransformComponent>().position();

		auto& grid = _marchingCubeComponent->grid();

		#pragma omp parallel for
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
							_marchingCubeComponent->changeGrid(x, y, z, i, pos);
						}
					}
				}
			}
		}
		std::uint64_t end = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		//std::cout << "Metaball: " << end-start << std::endl;
		
	}

	~MetaballComponent() override
	{
	}

private:
	MarchingCubeComponent* _marchingCubeComponent;
	float _radius;
};
