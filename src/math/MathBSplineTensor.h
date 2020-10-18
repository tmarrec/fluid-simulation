#pragma once

#include <cstdint>
#include <vector>
#include <tuple>

#include "../utils.h"
#include "MathBSpline.h"

class MathBSplineTensor
{
public:
	MathBSplineTensor(std::uint8_t __order, std::vector<std::vector<glm::vec3>> __controls)
	: _order { __order }
	, _controls { __controls }
	, _rangeV { (double)__order-1, (double)__controls.size() }
	{
		// Leading BSplines generation
		for (auto bs : __controls)
		{
			_leadingMathBSplines.emplace_back(MathBSpline{__order, bs});
		}
	}
	
	MathBSplineTensor(std::uint8_t __order, std::vector<std::vector<glm::vec3>> __controls, std::vector<std::vector<float>> __knots)
	: _order { __order }
	, _controls { __controls }
	, _rangeV { (double)__order-1, (double)__controls.size() }
	{
		// Leading BSplines generation
		ASSERT(__controls.size() == __knots.size(), "__controls and __knots needs to have the same size");
		for (std::uint64_t i = 0; i < __controls.size(); ++i)
		{
			_leadingMathBSplines.emplace_back(MathBSpline{__order, _controls[i], __knots[i]});
		}
	}

	glm::vec3 eval(float __u, float __v) const
	{
		std::vector<glm::vec3> generatorControls;
		for (auto lbs : _leadingMathBSplines)
		{
			generatorControls.emplace_back(lbs.eval(__u));
		}
		return MathBSpline(_order, generatorControls).eval(__v);
	}

	Range rangeU() const
	{
		auto minStart = _leadingMathBSplines.front().range().start;
		auto minEnd = _leadingMathBSplines.front().range().end;
		for (auto e : _leadingMathBSplines)
		{
			if (e.range().start > minStart)
			{
				minStart = e.range().start;
			}
			if (e.range().end < minEnd)
			{
				minEnd = e.range().end;
			}
		}
		return {minStart, minEnd};
	}
	Range rangeV() const { return _rangeV; }
	~MathBSplineTensor() {};

private:
	std::uint8_t _order;
	std::vector<std::vector<glm::vec3>> _controls;
	std::vector<MathBSpline> _leadingMathBSplines;
	const Range _rangeV;
};

