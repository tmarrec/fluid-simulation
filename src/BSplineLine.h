#pragma once

#include <algorithm>
#include <iterator>
#include <memory>
#include <tuple>
#include <numeric>

#include "math/MathBSpline.h"
#include "shapes.h"

class BSplineLine
{
public:
	BSplineLine(std::uint8_t __order, std::vector<glm::vec3> __controls, bool __show_controls, float __delta)
	: _bspline { std::make_unique<MathBSpline>(MathBSpline{__order, __controls}) }
	, _delta { __delta }
	, _show_controls { __show_controls }
	, _controls { __controls }
	{}

	Shape shape() const
	{
		std::vector<GLfloat> vertices;
		std::vector<GLfloat> normals;
		std::vector<GLuint> indices;
		auto range = _bspline->range();
		GLuint ind = 0;
			
		if (_show_controls)
		{
			std::vector<float> flat;
			for (auto p : _controls)
			{
				flat.emplace_back(p.x);
				flat.emplace_back(p.y);
				flat.emplace_back(p.z);
			}	
			vertices = flat;
			normals.assign(flat.size(), 1);
			indices.resize(_controls.size());
			std::iota(std::begin(indices), std::end(indices), 0);
			ind = indices.size()-1;
		}

		for (float i = range.start; i < range.end; i += _delta)
		{
			auto point = _bspline->eval(i);
			vertices.emplace_back(point.x);
			vertices.emplace_back(point.y);
			vertices.emplace_back(point.z);
			normals.insert(normals.end(), {1,1,1});
			indices.emplace_back(ind++);
		}
		
		Shape shape;
		shape.vertices = vertices;
		shape.normals = normals;
		shape.indices = indices;
		return shape;
	}

	~BSplineLine() {}

private:
	std::unique_ptr<MathBSpline> _bspline;
	float _delta;
	bool _show_controls;
	std::vector<glm::vec3> _controls;
};
