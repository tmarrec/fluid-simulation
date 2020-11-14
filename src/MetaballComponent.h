#pragma once

#include "ECS.h"
#include "TransformComponent.h"
#include "MarchingCubeComponent.h"
#include "glm/gtx/string_cast.hpp"
#include <cstdint>
#include <functional>

class MetaballComponent : public Component
{
public:
	MetaballComponent(MarchingCubeComponent* __marchingCubeComponent, float __radius)
	: Component{}
	, _marchingCubeComponent { __marchingCubeComponent }
	, _radius { __radius }
	{}

	void update([[maybe_unused]] double __deltaTime) override
	{
		// Testings
		t += 0.01f;
		if (entity->getEntityID() == 14)
		{
			auto& transformComponent = entity->getComponent<TransformComponent>();
			transformComponent.setPosition({0.0f, cos(t)*2.5f, sin(t)*2.5f});
		}
		else if (entity->getEntityID() == 15) 
		{
			auto& transformComponent = entity->getComponent<TransformComponent>();
			transformComponent.setPosition({0.0f, cos(t)*2.5f, -sin(t)*2.5f});
		}
		else
		{
			auto& transformComponent = entity->getComponent<TransformComponent>();
			transformComponent.setPosition({0.0f, sin(-t)*2.5f, 0.0f});
		}
		// End Testings

		ASSERT(entity->hasComponent<TransformComponent>(), "entity should have a TransformComponent");
		auto pos = entity->getComponent<TransformComponent>().position();
		auto& grid = _marchingCubeComponent->grid();
		std::function<float(glm::vec3)> func = std::bind(&MetaballComponent::f, pos, _radius, std::placeholders::_1); 
		_marchingCubeComponent->addFunc(func);

		//#pragma omp parallel for
		for (std::uint64_t x = 0; x < grid.size(); ++x)
		{
			for (std::uint64_t y = 0; y < grid[x].size(); ++y)
			{
				for (std::uint64_t z = 0; z < grid[x][y].size(); ++z)
				{
					// loop cell points	
					for (std::uint8_t i = 0; i < 8; ++i)
					{
						float inside = func(grid[x][y][z].points[i]);
						if (inside >= 0.5f)
						{
							_marchingCubeComponent->changeGrid(x, y, z, i, inside);
						}
					}
				}
			}
		}
	}
	
	static inline float f(glm::vec3 __center, float __radius, glm::vec3 __pos)
	{
		float d = glm::dot(__center - __pos, __center - __pos);
		if (d == 0.0f) return 0.0f;
		return std::pow(__radius, 2)/d;
	}

private:
	MarchingCubeComponent* _marchingCubeComponent;
	float _radius;
	double t = 0.0f;
};
