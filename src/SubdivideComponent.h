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

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

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
	SubdivideComponent()
	: Component{}
	{}

	void init()
	{
		ASSERT(entity->hasComponent<DrawableComponent>(), "entity should have a DrawableComponent");
		auto& drawableComponent = entity->getComponent<DrawableComponent>();
		auto vertices = drawableComponent.vertices();
		auto indices = drawableComponent.indices();
		ASSERT(vertices->size()%3 == 0, "vertices size should be power of 3");
		ASSERT(indices->size()%3 == 0, "indices size should be power of 3");

		_readMesh(vertices, indices);

		std::uint64_t start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		// SUBDIVIDE
		_subdivide(2);
		std::uint64_t end = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		std::cout << "Time : " << end-start << std::endl;

		_writeMesh(drawableComponent);

	}

private:
	void _subdivide(const std::uint64_t iterations)
	{
		_mesh.add_property(_edgePoint);
		_mesh.add_property(_vertexPoint);

		for (std::uint64_t i = 0; i < iterations; ++i)
		{
			std::unordered_map<Vertex_, MyMesh::VertexHandle> points;
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

			ASSERT(n >= 3, "point valency should be greater or equal to 3");
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

	std::vector<Vertex_> _getSubdividedVertices(std::vector<Vertex_>& vs, std::unordered_map<Vertex_, MyMesh::VertexHandle>& points)
	{
		// Find all faces vertices with properties on edges/points
		for (const auto& f_it : _mesh.faces())
		{
			// Each edges -> odd points
			for (auto fe_it = _mesh.fe_begin(f_it); fe_it.is_valid(); ++fe_it)
			{
				auto point = _mesh.property(_edgePoint, *fe_it);
				Vertex_ v = {point};
				if (points.find(v) == points.end())
				{
					points.insert({{v, _mesh.add_vertex(point)}});
				}
				vs.emplace_back(v);
			}

			// Each vertex -> even points 
			for (auto fv_it = _mesh.fv_begin(f_it); fv_it.is_valid(); ++fv_it)
			{
				auto point = _mesh.property(_vertexPoint, *fv_it);
				Vertex_ v = {point};
				if (points.find(v) == points.end())
				{
					points.insert({{v, _mesh.add_vertex(point)}});
				}
				vs.emplace_back(v);
			}
		}
		return vs;
	}


	void _updateFaces(std::vector<Vertex_>& vs, std::unordered_map<Vertex_, MyMesh::VertexHandle>& points)
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
		for (std::uint64_t i = 0; i < vs.size(); i += 6)
		{
			const auto& A = points.find(vs[i+0])->second; //         D        //
			const auto& B = points.find(vs[i+1])->second; //        /\        //
			const auto& C = points.find(vs[i+2])->second; //       /  \       //
			const auto& D = points.find(vs[i+3])->second; //      /    \      //
			const auto& E = points.find(vs[i+4])->second; //   A /_____ \ B   //
			const auto& F = points.find(vs[i+5])->second; //    / \    / \    //
			_mesh.add_face({A,B,C});                      //   /   \  /   \   //
			_mesh.add_face({D,B,A});                      //  /_____\/_____\  //
			_mesh.add_face({E,C,B});                      // F      C       E //
			_mesh.add_face({F,A,C});
		}
		_mesh.garbage_collection();
	}


	void _readMesh(std::shared_ptr<std::vector<GLfloat>> __vertices, std::shared_ptr<std::vector<GLuint>> __indices)
	{
		/*
		MyMesh::VertexHandle vhandle[__vertices->size()];
		for (std::uint64_t i = 0; i < __vertices->size(); i += 3)
		{
			std::cout << i << std::endl;
			vhandle[i/3] = _mesh.add_vertex(MyMesh::Point(__vertices->at(i), __vertices->at(i+1), __vertices->at(i+2)));
		}
		std::vector<MyMesh::VertexHandle> t;
		for (std::uint64_t i = 0; i < __indices->size(); i += 3)
		{
			t.clear();
			t.emplace_back(vhandle[__indices->at(i)]);
			t.emplace_back(vhandle[__indices->at(i+1)]);
			t.emplace_back(vhandle[__indices->at(i+2)]);
			_mesh.add_face(t);
		}
		std::cout << "Vertices : " << _mesh.n_vertices() << std::endl;
		std::cout << "Faces : " << _mesh.n_faces() << std::endl;
		*/
		MyMesh::VertexHandle vhandle[8];
		vhandle[0] = _mesh.add_vertex(MyMesh::Point(-1, -1,1));
		vhandle[1] = _mesh.add_vertex(MyMesh::Point( 1, -1,1));
		vhandle[2] = _mesh.add_vertex(MyMesh::Point( 1,1,1));
		vhandle[3] = _mesh.add_vertex(MyMesh::Point(-1,1,1));
		vhandle[4] = _mesh.add_vertex(MyMesh::Point(-1, -1, -1));
		vhandle[5] = _mesh.add_vertex(MyMesh::Point( 1, -1, -1));
		vhandle[6] = _mesh.add_vertex(MyMesh::Point( 1,1, -1));
		vhandle[7] = _mesh.add_vertex(MyMesh::Point(-1,1, -1));
		// generate (quadrilateral) faces
		std::vector<MyMesh::VertexHandle>face_vhandles;
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[0]);
		face_vhandles.push_back(vhandle[1]);
		face_vhandles.push_back(vhandle[2]);
		face_vhandles.push_back(vhandle[3]);
		_mesh.add_face(face_vhandles);

		face_vhandles.clear();
		face_vhandles.push_back(vhandle[7]);
		face_vhandles.push_back(vhandle[6]);
		face_vhandles.push_back(vhandle[5]);
		face_vhandles.push_back(vhandle[4]);
		_mesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[1]);
		face_vhandles.push_back(vhandle[0]);
		face_vhandles.push_back(vhandle[4]);
		face_vhandles.push_back(vhandle[5]);
		_mesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[2]);
		face_vhandles.push_back(vhandle[1]);
		face_vhandles.push_back(vhandle[5]);
		face_vhandles.push_back(vhandle[6]);
		_mesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[3]);
		face_vhandles.push_back(vhandle[2]);
		face_vhandles.push_back(vhandle[6]);
		face_vhandles.push_back(vhandle[7]);
		_mesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[0]);
		face_vhandles.push_back(vhandle[3]);
		face_vhandles.push_back(vhandle[7]);
		face_vhandles.push_back(vhandle[4]);
		_mesh.add_face(face_vhandles);
	}

	void _writeMesh(DrawableComponent& drawableComponent)
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
					newVertices.emplace_back(vTemp[0]);
					newVertices.emplace_back(vTemp[1]);
					newVertices.emplace_back(vTemp[2]);
					newNormals.emplace_back(nTemp[0]);
					newNormals.emplace_back(nTemp[1]);
					newNormals.emplace_back(nTemp[2]);
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

		drawableComponent.setVertices(newVertices);
		drawableComponent.setNormals(newNormals);
		drawableComponent.setIndices(newIndices);
		drawableComponent.updateGeometry();
	}

	MyMesh _mesh;
	OpenMesh::EPropHandleT<MyMesh::Point> _edgePoint;
	OpenMesh::VPropHandleT<MyMesh::Point> _vertexPoint;
};
