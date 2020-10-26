#pragma once

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Subdivider/Uniform/LoopT.hh>
#include <array>
#include <cstdint>
#include <unordered_map>

#include "ECS.h"
#include "DrawableComponent.h"
#include "OpenMesh/Core/Geometry/Vector11T.hh"
#include "OpenMesh/Core/Geometry/VectorT.hh"

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

struct Vertex_
{
	OpenMesh::VectorT<float, 3> first;
	bool operator==(const Vertex_ &other) const
	{ 
		return (first == other.first);
	}
};

namespace std
{
	template <>
	struct hash<Vertex_>
	{
		size_t operator()( const Vertex_& k ) const
		{
			// Compute individual hash values for first, second and third
			// http://stackoverflow.com/a/1646913/126995
			size_t res = 17;
			res = res * 31 + hash<float>()( k.first[0] );
			res = res * 31 + hash<float>()( k.first[1] );
			res = res * 31 + hash<float>()( k.first[2] );
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

		// SUBDIVIDE
		OpenMesh::Subdivider::Uniform::LoopT<MyMesh> l;
		l.attach(_mesh);
		l(1);
		l.detach();

		writeMesh(drawableComponent);
	}

private:
	void readMesh(std::shared_ptr<std::vector<GLfloat>> __vertices, std::shared_ptr<std::vector<GLuint>> __indices)
	{
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
	}

	void writeMesh(DrawableComponent& drawableComponent)
	{
		std::vector<GLfloat> newVertices;
		std::vector<GLuint> newIndices;
		std::unordered_map<Vertex_, std::uint64_t> points;
		std::uint64_t vi = 0;
		for (MyMesh::FaceIter f_it = _mesh.faces_sbegin(); f_it != _mesh.faces_end(); ++f_it)
		{
			GLuint ind[3];
			GLuint i = 0;
			for (MyMesh::ConstFaceHalfedgeIter fv_it = _mesh.cfh_iter(*f_it); fv_it.is_valid(); ++fv_it)
			{
				ASSERT(i < 3, "all faces needs to be triangles");
				auto p = _mesh.point(_mesh.to_vertex_handle(*fv_it));
				Vertex_ v = {p};
				auto got = points.find(v);
				if (got == points.end())
				{
					points.insert({{v, vi}});
					newVertices.emplace_back(p[0]);
					newVertices.emplace_back(p[1]);
					newVertices.emplace_back(p[2]);
					ind[i] = vi;
					vi++;
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
		drawableComponent.setVertices(newVertices);
		drawableComponent.setIndices(newIndices);
	}

	MyMesh _mesh;
};
