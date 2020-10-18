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

	glm::vec3 eval(float __u, float __v) const
	{
		std::vector<glm::vec3> generatorControls;
		for (auto lbs : _leadingMathBSplines)
		{
			generatorControls.emplace_back(lbs.eval(__u));
		}
		return MathBSpline(_order, generatorControls).eval(__v);
	}

	Range rangeU() const { return _leadingMathBSplines.front().range(); }
	Range rangeV() const { return _rangeV; }
	//std::uint64_t nbGenerato() const { return 
	~MathBSplineTensor() {};

private:
	std::uint8_t _order;
	std::vector<std::vector<glm::vec3>> _controls;
	std::vector<MathBSpline> _leadingMathBSplines;
	const Range _rangeV;
};

