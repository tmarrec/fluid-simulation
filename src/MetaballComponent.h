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
		float cellSize = 0.5f;
		// Grid initialization loop in middle of cell
		for (float x = -_radius*2; x < _radius*2; x += cellSize)
		{
			std::vector<std::vector<Cell>> grid2D;
			for (float y = -_radius*2; y < _radius*2; y += cellSize)
			{
				std::vector<Cell> grid1D;
				for (float z = -_radius*2; z < _radius*2; z += cellSize)
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
		std::vector<GLfloat> vertices;
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
					}
					std::uint64_t cubeIndex = 0;
					if (_grid[x][y][z].val[0] < isoLevel) cubeIndex |= 1;
					if (_grid[x][y][z].val[1] < isoLevel) cubeIndex |= 2;
					if (_grid[x][y][z].val[2] < isoLevel) cubeIndex |= 4;
					if (_grid[x][y][z].val[3] < isoLevel) cubeIndex |= 8;
					if (_grid[x][y][z].val[4] < isoLevel) cubeIndex |= 16;
					if (_grid[x][y][z].val[5] < isoLevel) cubeIndex |= 32;
					if (_grid[x][y][z].val[6] < isoLevel) cubeIndex |= 64;
					if (_grid[x][y][z].val[7] < isoLevel) cubeIndex |= 128;

					std::array<glm::vec3, 12> vertList;
					if (edgeTable[cubeIndex] & 1) vertList[0] = VertexInterp(isoLevel, _grid[x][y][z].points[0], _grid[x][y][z].points[1], _grid[x][y][z].val[0], _grid[x][y][z].val[1]);
					if (edgeTable[cubeIndex] & 2) vertList[1] = VertexInterp(isoLevel, _grid[x][y][z].points[1], _grid[x][y][z].points[2], _grid[x][y][z].val[1], _grid[x][y][z].val[2]);
					if (edgeTable[cubeIndex] & 4) vertList[2] = VertexInterp(isoLevel, _grid[x][y][z].points[2], _grid[x][y][z].points[3], _grid[x][y][z].val[2], _grid[x][y][z].val[3]);
					if (edgeTable[cubeIndex] & 8) vertList[3] = VertexInterp(isoLevel, _grid[x][y][z].points[3], _grid[x][y][z].points[0], _grid[x][y][z].val[3], _grid[x][y][z].val[0]);
					if (edgeTable[cubeIndex] & 16) vertList[4] = VertexInterp(isoLevel, _grid[x][y][z].points[4], _grid[x][y][z].points[5], _grid[x][y][z].val[4], _grid[x][y][z].val[5]);
					if (edgeTable[cubeIndex] & 32) vertList[5] = VertexInterp(isoLevel, _grid[x][y][z].points[5], _grid[x][y][z].points[6], _grid[x][y][z].val[5], _grid[x][y][z].val[6]);
					if (edgeTable[cubeIndex] & 64) vertList[6] = VertexInterp(isoLevel, _grid[x][y][z].points[6], _grid[x][y][z].points[7], _grid[x][y][z].val[6], _grid[x][y][z].val[7]);
					if (edgeTable[cubeIndex] & 128) vertList[7] = VertexInterp(isoLevel, _grid[x][y][z].points[7], _grid[x][y][z].points[4], _grid[x][y][z].val[7], _grid[x][y][z].val[4]);
					if (edgeTable[cubeIndex] & 256) vertList[8] = VertexInterp(isoLevel, _grid[x][y][z].points[0], _grid[x][y][z].points[4], _grid[x][y][z].val[0], _grid[x][y][z].val[4]);
					if (edgeTable[cubeIndex] & 512) vertList[9] = VertexInterp(isoLevel, _grid[x][y][z].points[1], _grid[x][y][z].points[5], _grid[x][y][z].val[1], _grid[x][y][z].val[5]);
					if (edgeTable[cubeIndex] & 1024) vertList[10] = VertexInterp(isoLevel, _grid[x][y][z].points[2], _grid[x][y][z].points[6], _grid[x][y][z].val[2], _grid[x][y][z].val[6]);
					if (edgeTable[cubeIndex] & 2048) vertList[11] = VertexInterp(isoLevel, _grid[x][y][z].points[3], _grid[x][y][z].points[7], _grid[x][y][z].val[3], _grid[x][y][z].val[7]);

					for (std::uint64_t i = 0; triTable[cubeIndex][i] != -1; i += 3)
					{
						glm::vec3 p1 = vertList[triTable[cubeIndex][i]];	
						glm::vec3 p2 = vertList[triTable[cubeIndex][i+1]];	
						glm::vec3 p3 = vertList[triTable[cubeIndex][i+2]];	
						vertices.emplace_back(p1.x);
						vertices.emplace_back(p1.y);
						vertices.emplace_back(p1.z);
						vertices.emplace_back(p2.x);
						vertices.emplace_back(p2.y);
						vertices.emplace_back(p2.z);
						vertices.emplace_back(p3.x);
						vertices.emplace_back(p3.y);
						vertices.emplace_back(p3.z);
					}
				}
			}
		}
		auto& drawable = entity->getComponent<DrawableComponent>();
		drawable.setVertices(vertices);
		ASSERT(vertices.size()%3==0, "vertices must be power of 3");
		std::vector<GLuint> indices;
		indices.resize(vertices.size()/3);
		std::iota(indices.begin(), indices.end(), 0);
		drawable.setIndices(indices);
		drawable.updateGeometry();
	}
	void draw() override
	{
	}

	void update([[maybe_unused]] double _deltaTime) override
	{
		// TODO ici modifier vertices
	}
	~MetaballComponent() override
	{
	}

private:
	float _radius;
	std::vector<std::vector<std::vector<Cell>>> _grid;
};
