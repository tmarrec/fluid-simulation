#pragma once

#include <GL/gl.h>
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <memory>
#include <tuple>
#include <numeric>

#include "glm/glm/geometric.hpp"
#include "math/MathBSplineTensor.h"
#include "shapes.h"

class BSplineTensor
{
public:
	BSplineTensor(std::uint8_t __order, std::vector<std::vector<glm::vec3>> __controls, bool __show_controls, float __delta)
	: _bsplineTensor { std::make_unique<MathBSplineTensor>(MathBSplineTensor{__order, __controls}) }
	, _delta { __delta }
	, _show_controls { __show_controls }
	, _controls { __controls }
	{}

	Shape shape() const
	{
		std::vector<GLfloat> vertices;
		std::vector<GLuint> indices;
		auto rangeU = _bsplineTensor->rangeU();
		auto rangeV = _bsplineTensor->rangeV();

		std::uint64_t nbPointsGenerator = 0;
		for (float u = rangeU.start; u < rangeU.end; u += _delta)
			nbPointsGenerator++;
		std::uint64_t nbGenerator = 0;
		for (float v = rangeV.start; v < rangeV.end; v += _delta)
			nbGenerator++;

		std::vector<glm::vec3> points;
		for (float u = rangeU.start; u < rangeU.end; u += _delta)
		{
			for (float v = rangeV.start; v < rangeV.end; v += _delta)
			{
				auto point = _bsplineTensor->eval(u, v);
				points.emplace_back(point);
				vertices.emplace_back(point.x);
				vertices.emplace_back(point.y);
				vertices.emplace_back(point.z);
			}
		}
		std::vector<GLfloat> normals(vertices.size());
		std::vector<std::vector<glm::vec3>> normalsPerPoint(points.size());

		std::cout << nbPointsGenerator << " " << nbGenerator << std::endl;
		for (std::uint64_t i = 0; i < nbPointsGenerator-1; ++i)
		{
			for (std::uint64_t j = 0; j < nbGenerator-1; ++j)
			{
				GLuint p = j+i*nbGenerator;	
				glm::vec3 norm;

				// First face triangle
				indices.emplace_back(p);
				indices.emplace_back(p+nbGenerator);
				indices.emplace_back(p+1);

				norm = glm::normalize(glm::cross(points[p+nbGenerator]-points[p], points[p+1]-points[p]));
				normalsPerPoint[p].emplace_back(norm);
				normalsPerPoint[p+nbGenerator].emplace_back(norm);
				normalsPerPoint[p+1].emplace_back(norm);

				std::cout << p << " " << p+nbGenerator << " " << p+1 << std::endl;

				// Second face triangle
				indices.emplace_back(p+nbGenerator);
				indices.emplace_back(p+nbGenerator+1);
				indices.emplace_back(p+1);

				norm = glm::normalize(glm::cross(points[p+nbGenerator+1]-points[p+nbGenerator], points[p+1]-points[p+nbGenerator]));
				normalsPerPoint[p+nbGenerator].emplace_back(norm);
				normalsPerPoint[p+nbGenerator+1].emplace_back(norm);
				normalsPerPoint[p+1].emplace_back(norm);
				std::cout << p+nbGenerator << " " << p+nbGenerator+1 << " " << p+1 << std::endl << std::endl;
			}
		}

		for (std::uint64_t i = 0; i < normalsPerPoint.size(); ++i)
		{
			glm::vec3 average = {0.0f, 0.0f, 0.0f};
			for (auto a : normalsPerPoint[i])
			{
				average += a;
			}
			average /= normalsPerPoint[i].size();
			normals[3*i] = average.x;
			normals[3*i+1] = average.y;
			normals[3*i+2] = average.z;
		}
		std::cout << std::endl;

		Shape shape;
		shape.vertices = vertices;
		shape.normals = normals;
		shape.indices = indices;
		return shape;
	}

	~BSplineTensor() {}

private:
	std::unique_ptr<MathBSplineTensor> _bsplineTensor;
	float _delta;
	bool _show_controls;
	std::vector<std::vector<glm::vec3>> _controls;
};
