#pragma once

#include <GL/gl.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Subdivider/Uniform/LoopT.hh>
#include <array>
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
		/*
		OpenMesh::Subdivider::Uniform::LoopT<MyMesh> l;
		l.attach(_mesh);
		l(2);
		l.detach();
		*/
		_subdivide();
		std::uint64_t end = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		std::cout << "Time : " << end-start << std::endl;

		_writeMesh(drawableComponent);

	}

private:
	void _subdivide()
	{
		OpenMesh::EPropHandleT<MyMesh::Point> edgePoint;
		OpenMesh::VPropHandleT<MyMesh::Point> vertexPoint;
		_mesh.add_property(edgePoint);
		_mesh.add_property(vertexPoint);

		// Odd Vertices
		for (auto e_it = _mesh.edges_begin(); e_it != _mesh.edges_end(); ++e_it)
		{
			auto edge = _mesh.edge(*e_it);
			auto v1 = _mesh.point(_mesh.to_vertex_handle(_mesh.halfedge_handle(*e_it, 0)));
			auto v2 = _mesh.point(_mesh.from_vertex_handle(_mesh.halfedge_handle(*e_it, 0)));
			//std::cout << v1 << " -> " << v2 << " -> " << (v2+v1)/2 << std::endl;
			_mesh.property(edgePoint, *e_it) = (v2+v1)/2;
		}
		// Even Vertices
		for (auto v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
		{
			auto point = _mesh.point(*v_it);
			/*
			auto n = _mesh.valence(*v_it);	
			ASSERT(n >= 3, "point valance should be greater or equal to 3");
			float B;
			if (n > 3)
			{
				B = 3/(8*(float)n);
			}
			else
			{
				B = 3/16;
			}
			std::cout << point << std::endl;	
			std::cout << n << " " << B << std::endl;
			std::cout << std::endl;
			*/


			//
			_mesh.property(vertexPoint, *v_it) = point;
		}

		
		std::vector<MyMesh::VertexHandle> vhandle;
		std::vector<std::vector<MyMesh::VertexHandle>> fHandle;

		// New faces creation
		for (auto f_it = _mesh.faces_begin(); f_it != _mesh.faces_end(); ++f_it)
		{
			// Each edges
			std::vector<MyMesh::VertexHandle> veHandle;
			for (auto fe_it = _mesh.fe_begin(*f_it); fe_it.is_valid(); ++fe_it)
			{
				auto point = _mesh.property(edgePoint, *fe_it);
				veHandle.emplace_back(_mesh.add_vertex(point));
			}

			// Each vertex
			std::vector<MyMesh::VertexHandle> vvHandle;
			for (auto fv_it = _mesh.fv_begin(*f_it); fv_it.is_valid(); ++fv_it)
			{
				auto point = _mesh.property(vertexPoint, *fv_it);
				vvHandle.emplace_back(_mesh.add_vertex(point));
			}
		
			fHandle.emplace_back(veHandle);

			std::vector<MyMesh::VertexHandle> test;
			test.emplace_back(vvHandle[0]);
			test.emplace_back(veHandle[1]);
			test.emplace_back(veHandle[0]);
			fHandle.emplace_back(test);
			test.clear();
			test.emplace_back(vvHandle[1]);
			test.emplace_back(veHandle[2]);
			test.emplace_back(veHandle[1]);
			fHandle.emplace_back(test);
			test.clear();
			test.emplace_back(vvHandle[2]);
			test.emplace_back(veHandle[0]);
			test.emplace_back(veHandle[2]);
			fHandle.emplace_back(test);
			test.clear();

		}
		_mesh.request_face_status();
		_mesh.request_edge_status();
		_mesh.request_vertex_status();
		for (auto f_it = _mesh.faces_begin(); f_it != _mesh.faces_end(); ++f_it)
		{
			_mesh.delete_face(*f_it, true);
		}
		// Create all faces
		for (const auto& fh : fHandle)
		{
			_mesh.add_face(fh);
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
};
