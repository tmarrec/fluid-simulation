#pragma once

#include <GL/gl.h>
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <memory>
#include <tuple>
#include <numeric>

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
		std::vector<GLfloat> normals;
		std::vector<GLuint> indices;
		auto rangeU = _bsplineTensor->rangeU();
		auto rangeV = _bsplineTensor->rangeV();

		std::uint64_t nbPointsGenerator = 0;
		for (float u = rangeU.start; u < rangeU.end; u += _delta)
			nbPointsGenerator++;
		std::uint64_t nbGenerator = 0;
		for (float v = rangeV.start; v < rangeV.end; v += _delta)
			nbGenerator++;

		for (float u = rangeU.start; u < rangeU.end; u += _delta)
		{
			for (float v = rangeV.start; v < rangeV.end; v += _delta)
			{
				auto point = _bsplineTensor->eval(u, v);
				vertices.emplace_back(point.x);
				vertices.emplace_back(point.y);
				vertices.emplace_back(point.z);
				normals.emplace_back(1);
				normals.emplace_back(1);
				normals.emplace_back(1);
			}
		}
		std::cout << nbPointsGenerator << " " << nbGenerator << std::endl;
		for (std::uint64_t i = 0; i < nbPointsGenerator-1; ++i)
		{
			for (std::uint64_t j = 0; j < nbGenerator-1; ++j)
			{
				std::cout << i << " " << j << std::endl;
				GLuint p = j+i*nbGenerator;	

				// First face triangle
				indices.emplace_back(p);
				indices.emplace_back(p+nbGenerator);
				indices.emplace_back(p+1);
				std::cout << p << " " << p+nbGenerator << " " << p+1 << std::endl;

				// Second face triangle
				indices.emplace_back(p+nbGenerator);
				indices.emplace_back(p+nbGenerator+1);
				indices.emplace_back(p+1);
				std::cout << p+nbGenerator << " " << p+nbGenerator+1 << " " << p+1 << std::endl << std::endl;
			}
		}

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
