#pragma once

#include <GL/gl.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Subdivider/Uniform/LoopT.hh>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <chrono>

#include "ECS.h"
#include "DrawableComponent.h"
#include "OpenMesh/Core/Geometry/Vector11T.hh"
#include "OpenMesh/Core/Geometry/VectorT.hh"
#include "OpenMesh/Core/Mesh/PolyConnectivity.hh"
#include "OpenMesh/Core/Utils/Property.hh"

typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

// Vertex struct with hash function used for hash map
struct Vertex_
{
	OpenMesh::VectorT<float, 3> point;
	bool operator==(const Vertex_ &other) const
	{ 
		return (point == other.point);
	}
};

namespace std
{
	template <>
	struct hash<Vertex_>
	{
		size_t operator()(const Vertex_& k) const
		{
			size_t res = 17;
			res = res * 31 + hash<float>()(k.point[0]);
			res = res * 31 + hash<float>()(k.point[1]);
			res = res * 31 + hash<float>()(k.point[2]);
			return res;
		}
	};
}

class SubdivideComponent : public Component
{
public:
	explicit SubdivideComponent(float __n)
	: Component{}
	, _n { __n }
	{}

	void init() override
	{
		ASSERT(entity->hasComponent<DrawableComponent>(), "entity should have a DrawableComponent");
		auto& drawableComponent = entity->getComponent<DrawableComponent>();
		auto vertices = drawableComponent.vertices();
		auto indices = drawableComponent.indices();
		ASSERT(vertices->size()%3 == 0, "vertices size should be power of 3");
		ASSERT(indices->size()%3 == 0, "indices size should be power of 3");

		_readMesh(vertices, indices);
		_subdivide(_n);
		_writeMesh(drawableComponent);
	}

private:
	float _n;

	void _subdivide(const std::uint64_t __iterations)
	{
		_mesh.add_property(_edgePoint);
		_mesh.add_property(_vertexPoint);

		for (std::uint64_t i = 0; i < __iterations; ++i)
		{
			std::unordered_map<Vertex_, Mesh::VertexHandle> points;
			std::vector<Vertex_> vs;
			_smoothPoints();
			_getSubdividedVertices(vs, points);
			_updateFaces(vs, points);
		}
	}

	void _smoothPoints() 
	{
		// Odd Vertices
		for (const auto& e_it : _mesh.edges())
		{
			auto he = _mesh.halfedge_handle(e_it, 0);
			auto heOp = _mesh.next_halfedge_handle(he);
			auto heNext = _mesh.next_halfedge_handle(heOp);
			auto a = _mesh.point(_mesh.to_vertex_handle(he));
			auto b = _mesh.point(_mesh.to_vertex_handle(heNext));
			OpenMesh::VectorT<float, 3> point;

			if (_mesh.is_boundary(e_it))
			{
				point = (a+b)/2.0f;
			}
			else
			{
				auto c = _mesh.point(_mesh.to_vertex_handle(heOp));
				auto heOp2 = _mesh.next_halfedge_handle(_mesh.halfedge_handle(e_it, 1));
				auto d = _mesh.point(_mesh.to_vertex_handle(heOp2));
				point = 3.0f/8.0f*(a+b)+1.0f/8.0f*(c+d);
			}
			_mesh.property(_edgePoint, e_it) = point;
		}
		// Even Vertices
		for (const auto& v_it : _mesh.vertices())
		{
			auto point = _mesh.point(v_it);
			auto n = _mesh.valence(v_it);	
			float B;

			ASSERT(n >= 3, "point valence should be greater or equal to 3");
			if (n > 3)
			{
				B = (1.0f/(float)n)*((5.0f/8.0f)-std::pow((3.0f/8.0f)+((1.0f/4.0f)*cos((2.0f*M_PIf64)/(float)n)), 2));
			}
			else
			{
				B = 3.0f/16.0f;
			}
			OpenMesh::VectorT<float,3> sum = {0, 0, 0}; 
			for (auto vv_it = _mesh.vv_iter(v_it); vv_it.is_valid(); ++vv_it)
			{
				sum += _mesh.point(*vv_it);
			}
			point = point*(1-n*B)+(sum*B);
			_mesh.property(_vertexPoint, v_it) = point;
		}
	}

	std::vector<Vertex_> _getSubdividedVertices(std::vector<Vertex_>& __vs, std::unordered_map<Vertex_, Mesh::VertexHandle>& __points)
	{
		// Find all faces vertices with properties on edges/points
		for (const auto& f_it : _mesh.faces())
		{
			// Each edges -> add odd points
			for (auto fe_it = _mesh.fe_begin(f_it); fe_it.is_valid(); ++fe_it)
			{
				auto point = _mesh.property(_edgePoint, *fe_it);
				Vertex_ v = {point};
				if (__points.find(v) == __points.end())
				{
					__points.insert({{v, _mesh.add_vertex(point)}});
				}
				__vs.emplace_back(v);
			}

			// Each vertex -> add even points 
			for (auto fv_it = _mesh.fv_begin(f_it); fv_it.is_valid(); ++fv_it)
			{
				auto point = _mesh.property(_vertexPoint, *fv_it);
				Vertex_ v = {point};
				if (__points.find(v) == __points.end())
				{
					__points.insert({{v, _mesh.add_vertex(point)}});
				}
				__vs.emplace_back(v);
			}
		}
		return __vs;
	}


	void _updateFaces(std::vector<Vertex_>& __vs, std::unordered_map<Vertex_, Mesh::VertexHandle>& __points)
	{
		_mesh.request_face_status();
		_mesh.request_edge_status();
		_mesh.request_vertex_status();
		// Delete all faces
		for (const auto& f_it : _mesh.faces())
		{
			_mesh.delete_face(f_it, true);
		}
		// Create all new faces
		for (std::uint64_t i = 0; i < __vs.size(); i += 6)
		{
			const auto& A = __points.find(__vs[i+0])->second; //         D        //
			const auto& B = __points.find(__vs[i+1])->second; //        /\        //
			const auto& C = __points.find(__vs[i+2])->second; //       /  \       //
			const auto& D = __points.find(__vs[i+3])->second; //      /    \      //
			const auto& E = __points.find(__vs[i+4])->second; //   A /_____ \ B   //
			const auto& F = __points.find(__vs[i+5])->second; //    / \    / \    //
			_mesh.add_face({A,B,C});                          //   /   \  /   \   //
			_mesh.add_face({D,B,A});                          //  /_____\/_____\  //
			_mesh.add_face({E,C,B});                          // F      C       E //
			_mesh.add_face({F,A,C});
		}
		_mesh.garbage_collection();
	}


	void _readMesh(std::shared_ptr<std::vector<GLfloat>> __vertices, std::shared_ptr<std::vector<GLuint>> __indices)
	{
		std::unordered_map<Vertex_, Mesh::VertexHandle> points;
		std::vector<Mesh::VertexHandle> e;

		Mesh::VertexHandle vhandle[__vertices->size()];
		for (std::uint64_t i = 0; i < __vertices->size(); i += 3)
		{
			Mesh::Point p (__vertices->at(i), __vertices->at(i+1), __vertices->at(i+2));
			Vertex_ v {p};
			auto got = points.find(v);
			if (got == points.end())
			{
				auto k = _mesh.add_vertex(Mesh::Point(__vertices->at(i), __vertices->at(i+1), __vertices->at(i+2)));
				points.insert({{v, k}});
				e.emplace_back(k);
			}
			else
			{
				e.emplace_back(got->second);
			}
		}
		std::vector<Mesh::VertexHandle> t {3};
		for (std::uint64_t i = 0; i < __indices->size(); i += 3)
		{
			t[0] = e.at(__indices->at(i));
			t[1] = e.at(__indices->at(i+1));
			t[2] = e.at(__indices->at(i+2));
			_mesh.add_face(t);
		}
	}

	void _writeMesh(DrawableComponent& __drawableComponent)
	{
		std::vector<GLfloat> newVertices;
		std::vector<GLfloat> newNormals;
		std::vector<GLuint> newIndices;
		std::unordered_map<Vertex_, std::uint64_t> points;
		std::uint64_t vertexInd = 0;
		_mesh.request_vertex_normals();
		_mesh.request_halfedge_normals();
		_mesh.request_face_normals();
		_mesh.update_normals();
		for (auto f_it = _mesh.faces_sbegin(); f_it != _mesh.faces_end(); ++f_it)
		{
			std::vector<std::uint64_t> ind {3};
			std::uint64_t i = 0;
			for (auto fv_it = _mesh.cfh_iter(*f_it); fv_it.is_valid(); ++fv_it)
			{
				ASSERT(i < 3, "all faces needs to be triangles");
				auto vTemp = _mesh.point(_mesh.to_vertex_handle(*fv_it));
				auto nTemp = _mesh.normal(_mesh.to_vertex_handle(*fv_it));
				Vertex_ v = {vTemp};
				auto got = points.find(v);
				if (got == points.end())
				{
					points.insert({{v, vertexInd}});
					std::uint64_t oldSize = newVertices.size();
					newVertices.resize(oldSize+3);
					memcpy(&newVertices[oldSize], &vTemp, sizeof(GLfloat)*3);
					newNormals.resize(oldSize+3);
					memcpy(&newNormals[oldSize], &nTemp, sizeof(GLfloat)*3);
					ind[i] = vertexInd++;
				}
				else
				{
					ind[i] = got->second;
				}
				i++;
			}
			newIndices.emplace_back(ind[0]);
			newIndices.emplace_back(ind[1]);
			newIndices.emplace_back(ind[2]);
		}
		_mesh.release_face_normals();
		_mesh.release_halfedge_normals();
		_mesh.release_vertex_normals();

		__drawableComponent.setVertices(newVertices);
		__drawableComponent.setNormals(newNormals);
		__drawableComponent.setIndices(newIndices);
		__drawableComponent.updateGeometry();
	}

	Mesh _mesh;
	OpenMesh::EPropHandleT<Mesh::Point> _edgePoint;
	OpenMesh::VPropHandleT<Mesh::Point> _vertexPoint;
};
