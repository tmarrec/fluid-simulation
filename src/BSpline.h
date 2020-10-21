#pragma once

#include <memory>
#include <numeric>

#include "shapes.h"

// Math representation of a BSpline
class BSpline
{
public:
	// BSpline without custom knots -> uniform knots
	BSpline(std::uint8_t __order, std::vector<glm::vec3> __controls)
	: _order { __order }
	, _knots { _uniform_vector(__controls.size()+__order+1) }
	, _controls { __controls }
	, _range { _knots[__order-1], _knots[__controls.size()] }
	{}
	// BSpline with custom knots
	BSpline(std::uint8_t __order, std::vector<glm::vec3> __controls, std::vector<float> __knots)
	: _order { __order }
	, _knots { __knots }
	, _controls { __controls }
	, _range { _knots[__order-1], _knots[__controls.size()] }
	{}

	// Evaluate the point on the BSpline for __u
	glm::vec3 eval(float __u) const
	{
		// Blooming
		ASSERT(__u >= _range.start && __u <= _range.end, "__u is not in the range of the BSpline");
		std::uint8_t dec = 0;
		std::uint8_t i = _order;
		std::uint8_t k = _order;
		std::uint8_t m = _order-1;
		std::vector<glm::vec3> P_temps;
		while (__u > _knots[i])
		{
			i++;
			dec++;
		}
		P_temps.insert(P_temps.end(), _controls.begin()+dec, _controls.begin()+dec+_order); 
		for (std::uint8_t l = 0; l < m; ++l)
		{
			for (std::uint8_t j = 0; j < k-1; ++j)
			{
				P_temps[j] = ((_knots[dec+k+j]-__u)/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j]+((__u-_knots[dec+1+j])/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j+1];
			}
			dec++;
			k--;
		}
		return P_temps[0];
	}

	Range range() const { return _range; }
	~BSpline() {};

private:
	// Uniform vector generation starting from zero
	std::vector<float> _uniform_vector(std::uint8_t __size) const
	{
		std::vector<float> v(__size);
		std::iota(v.begin(), v.end(), 0);
		return v;
	}

	std::uint8_t _order;
	std::vector<float> _knots;
	std::vector<glm::vec3> _controls;
	Range _range;
};

// Shape representation of a BSpline
class BSplineShape
{
public:
	// Unspecified knots -> uniform knots
	BSplineShape(std::uint8_t __order, std::vector<glm::vec3> __controls, bool __show_controls, float __delta)
	: _bspline { std::make_unique<BSpline>(BSpline{__order, __controls}) }
	, _delta { __delta }
	, _show_controls { __show_controls }
	, _controls { __controls }
	{}
	// Specified knots
	BSplineShape(std::uint8_t __order, std::vector<glm::vec3> __controls, bool __show_controls, float __delta, std::vector<float> __knots)
	: _bspline { std::make_unique<BSpline>(BSpline{__order, __controls, __knots}) }
	, _delta { __delta }
	, _show_controls { __show_controls }
	, _controls { __controls }
	{}

	// Generate the vertices, normals and indices of the BSpline shape with
	// the math representation of the BSpline
	Shape shape() const
	{
		std::vector<GLfloat> vertices;
		std::vector<GLfloat> normals;
		std::vector<GLuint> indices;
		auto range = _bspline->range();
		// Controls attributes
		if (_show_controls)
		{
			for (const auto& p : _controls)
			{
				vertices.insert(vertices.end(), {p.x, p.y, p.z}); 
			}	
		}
		// BSpline attributes
		for (float i = range.end; i > range.start; i -= _delta)
		{
			auto p = _bspline->eval(i);
			vertices.insert(vertices.end(), {p.x, p.y, p.z});
		}
		ASSERT(vertices.size()%3 == 0, "vertices should be power of 3");
		normals.assign(vertices.size(), 1);
		indices.resize((vertices.size()/3)-1);
		std::iota(indices.begin(), indices.end(), 0);
		return {vertices, normals, indices};
	}

	~BSplineShape() {}

private:
	std::unique_ptr<BSpline> _bspline;
	float _delta;
	bool _show_controls;
	std::vector<glm::vec3> _controls;
};
