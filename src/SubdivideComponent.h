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

		readMesh(vertices, indices);

		std::uint64_t start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		// SUBDIVIDE
		OpenMesh::Subdivider::Uniform::LoopT<MyMesh> l;
		l.attach(_mesh);
		l(2);
		l.detach();
		std::uint64_t end = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		std::cout << "Time : " << end-start << std::endl;

		writeMesh(drawableComponent);

	}

private:
	void readMesh(std::shared_ptr<std::vector<GLfloat>> __vertices, std::shared_ptr<std::vector<GLuint>> __indices)
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

	void writeMesh(DrawableComponent& drawableComponent)
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
		for (MyMesh::FaceIter f_it = _mesh.faces_sbegin(); f_it != _mesh.faces_end(); ++f_it)
		{
			std::vector<std::uint64_t> ind {3};
			std::uint64_t i = 0;
			for (MyMesh::ConstFaceHalfedgeIter fv_it = _mesh.cfh_iter(*f_it); fv_it.is_valid(); ++fv_it)
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
