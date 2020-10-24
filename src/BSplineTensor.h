#pragma once

#include "BSpline.h"

// Math representation of a BSplien Tensor
class BSplineTensor
{
public:
	// BSplineTensor without custom knots -> uniform knots
	BSplineTensor(std::uint8_t __order, std::vector<std::vector<glm::vec3>> __controls)
	: _order { __order }
	, _controls { __controls }
	, _rangeV { (double)__order-1, (double)__controls.size() }
	{
		// Leading BSplines generation
		for (const auto& bs : __controls)
		{
			_leadingBSplines.emplace_back(BSpline{__order, bs});
		}
	}
	// BSplineTensor with custom knots	
	BSplineTensor(std::uint8_t __order, std::vector<std::vector<glm::vec3>> __controls, std::vector<std::vector<float>> __knots)
	: _order { __order }
	, _controls { __controls }
	, _rangeV { (double)__order-1, (double)__controls.size() }
	{
		// Leading BSplines generation
		ASSERT(__controls.size() == __knots.size(), "__controls and __knots needs to have the same size");
		for (std::uint64_t i = 0; i < __controls.size(); ++i)
		{
			_leadingBSplines.emplace_back(BSpline{__order, _controls[i], __knots[i]});
		}
	}

	// Evaluate the point on the BSplineTensor for __u and __v
	glm::vec3 eval(float __u, float __v) const
	{
		std::vector<glm::vec3> generatorControls;
		for (const auto& lbs : _leadingBSplines)
		{
			generatorControls.emplace_back(lbs.eval(__u));
		}
		return BSpline(_order, generatorControls).eval(__v);
	}

	// Get the narrowest U range
	Range rangeU() const
	{
		auto minStart = _leadingBSplines[0].range().start;
		auto minEnd = _leadingBSplines[0].range().end;
		for (const auto& e : _leadingBSplines)
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
	~BSplineTensor() {};

private:
	std::uint8_t _order;
	std::vector<std::vector<glm::vec3>> _controls;
	std::vector<BSpline> _leadingBSplines;
	const Range _rangeV;
};

// Shape representation of a BSpline Tensor
class BSplineTensorShape
{
public:
	// No custom knots -> uniform knots
	BSplineTensorShape(std::uint8_t __order, std::vector<std::vector<glm::vec3>> __controls, float __delta)
	: _bsplineTensor { std::make_unique<BSplineTensor>(BSplineTensor{__order, __controls}) }
	, _delta { __delta }
	, _controls { __controls }
	{}
	// Custom knots
	BSplineTensorShape(std::uint8_t __order, std::vector<std::vector<glm::vec3>> __controls, float __delta, std::vector<std::vector<float>> __knots)
	: _bsplineTensor { std::make_unique<BSplineTensor>(BSplineTensor{__order, __controls, __knots}) }
	, _delta { __delta }
	, _controls { __controls }
	{}

	// Generate the vertices, normals and indices of the BSplineTensor shape with
	// the math representation of the BSplineTensor
	Shape shape() const
	{
		std::vector<GLfloat> vertices;
		std::vector<GLuint> indices;
		std::uint64_t nbPointsGenerator = 0;
		std::uint64_t nbGenerator = 0;
		auto rangeU = _bsplineTensor->rangeU();
		auto rangeV = _bsplineTensor->rangeV();

		// Vertices generation with the math representation
		// of the BSplineTensor
		std::vector<glm::vec3> points;
		for (float u = rangeU.start; u < rangeU.end; u += _delta)
		{
			nbPointsGenerator++;
			for (float v = rangeV.start; v < rangeV.end; v += _delta)
			{
				auto point = _bsplineTensor->eval(u, v);
				points.emplace_back(point);
				vertices.insert(vertices.end(), {point.x, point.y, point.z});
				if (u == rangeU.start)
				{
					nbGenerator++;
				}
			}
		}
		// Normals for each triangles and indices computation
		std::vector<GLfloat> normals(vertices.size());
		std::vector<std::vector<glm::vec3>> normalsPerPoint(points.size());
		for (std::uint64_t i = 0; i < nbPointsGenerator-1; ++i)
		{
			for (std::uint64_t j = 0; j < nbGenerator-1; ++j)
			{
				GLuint p = j+i*nbGenerator;	
				// For the two triangles : p+nbGenerator+1 -- p+nbGenerator
				//                                |    2    /   1    |
				//                               p+1 --------------- p
				// Triangles face normal
				glm::vec3 normTriangle1 = glm::normalize(glm::cross(points[p+nbGenerator]-points[p], points[p+1]-points[p]));
				glm::vec3 normTriangle2 = glm::normalize(glm::cross(points[p+nbGenerator+1]-points[p+nbGenerator], points[p+1]-points[p+nbGenerator]));
				
				// Add those normals to each points adjacent to one or more of thoses triangles 
				normalsPerPoint[p].emplace_back(normTriangle1);
				normalsPerPoint[p+1].insert(normalsPerPoint[p+1].end(), {normTriangle1, normTriangle2});
				normalsPerPoint[p+nbGenerator].insert(normalsPerPoint[p+nbGenerator].end(), {normTriangle1, normTriangle2});
				normalsPerPoint[p+nbGenerator+1].emplace_back(normTriangle2);
				
				// Add both triangles indices
				indices.insert(indices.end(), {p, GLuint(p+nbGenerator), p+1, GLuint(p+nbGenerator), GLuint(p+nbGenerator+1), p+1});
			}
		}
		// Compute normals for each points with the average of their adjacent triangle faces
		for (std::uint64_t i = 0; i < normalsPerPoint.size(); ++i)
		{
			glm::vec3 average = {0.0f, 0.0f, 0.0f};
			for (const auto& a : normalsPerPoint[i])
			{
				average += a;
			}
			average /= normalsPerPoint[i].size();
			normals[3*i] = -average.x;
			normals[3*i+1] = -average.y;
			normals[3*i+2] = -average.z;
		}
		return {vertices, normals, indices};
	}

	~BSplineTensorShape() {}

private:
	std::unique_ptr<BSplineTensor> _bsplineTensor;
	float _delta;
	std::vector<std::vector<glm::vec3>> _controls;
};
