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
		// Grid initialization loop in middle of cell
		for (float x = 0; x < 10; x += cellSize)
		{
			std::vector<std::vector<Cell>> grid2D;
			for (float y = 0; y < 10; y += cellSize)
			{
				std::vector<Cell> grid1D;
				for (float z = 0; z < 10; z += cellSize)
				{
					Cell cell;
					float s = cellSize/2;
					cell.points[0] = {x-s,y+s,z-s}; 
					cell.points[1] = {x+s,y+s,z-s}; 
					cell.points[2] = {x+s,y-s,z-s}; 
					cell.points[3] = {x-s,y-s,z-s}; 
					cell.points[4] = {x-s,y+s,z+s}; 
					cell.points[5] = {x+s,y+s,z+s}; 
					cell.points[6] = {x+s,y-s,z+s}; 
					cell.points[7] = {x-s,y-s,z+s}; 
					for (std::uint64_t i = 0; i < 8; ++i)
					{
						cell.val[i] = 0;
					}
					grid1D.emplace_back(cell);
				}
				grid2D.emplace_back(grid1D);
			}
			_grid.emplace_back(grid2D);
		}
		
		// TODO do this at update
		// check inside grid
		auto pos = entity->getComponent<TransformComponent>().position();
		float isoLevel = 1.0f;
		for (std::uint64_t x = 0; x < _grid.size(); ++x)
		{
			for (std::uint64_t y = 0; y < _grid[x].size(); ++y)
			{
				for (std::uint64_t z = 0; z < _grid[x][y].size(); ++z)
				{
					// loop cell points	
					for (std::uint64_t i = 0; i < 8; ++i)
					{
						float inside = std::pow(_radius, 2)/(std::pow(_grid[x][y][z].points[i].x-pos.x, 2)+std::pow(_grid[x][y][z].points[i].y-pos.y, 2)+std::pow(_grid[x][z][z].points[i].z-pos.z, 2));
						if (inside >= 1.0f)
						{
							_grid[x][y][z].val[i] = isoLevel;
						}
						std::cout << _grid[x][y][z].val[i] << " ";
					}
					std::cout << std::endl;
					std::uint64_t cubeIndex = 0;
					if (_grid[x][y][z].val[0] < isoLevel) cubeIndex |= 1;
					if (_grid[x][y][z].val[1] < isoLevel) cubeIndex |= 2;
					if (_grid[x][y][z].val[2] < isoLevel) cubeIndex |= 4;
					if (_grid[x][y][z].val[3] < isoLevel) cubeIndex |= 8;
					if (_grid[x][y][z].val[4] < isoLevel) cubeIndex |= 16;
					if (_grid[x][y][z].val[5] < isoLevel) cubeIndex |= 32;
					if (_grid[x][y][z].val[6] < isoLevel) cubeIndex |= 64;
					if (_grid[x][y][z].val[7] < isoLevel) cubeIndex |= 128;
					std::cout << cubeIndex << std::endl << std::endl;

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
	std::vector<std::vector<std::vector<Cell>>> _grid;
};
