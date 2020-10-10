#pragma once

#include <algorithm>
#include <iterator>
#include <memory>
#include <tuple>
#include <numeric>

#include "math/BSpline.h"
#include "shapes.h"

class BSplineLine
{
public:
	BSplineLine(std::uint8_t __order, std::vector<glm::vec3> __controls, bool __show_controls, float __delta)
	: _bspline { std::make_unique<BSpline>(BSpline{__order, uniform_vector(__controls.size()+__order+1), __controls}) }
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
			
		if (_show_controls) {
			std::vector<float> flat;
			for (auto p : _controls) {
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

		for (float i = range.start; i < range.end; i += _delta) {
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
	std::unique_ptr<BSpline> _bspline;
	float _delta;
	bool _show_controls;
	std::vector<glm::vec3> _controls;

	std::vector<float> uniform_vector(std::uint8_t __size) const
	{
		std::vector<float> v;
		std::generate_n(std::back_inserter(v), __size,
		[]() -> float
		{
			static std::uint8_t i = 0;
			return i++;
		});
		return v;
	}
};
