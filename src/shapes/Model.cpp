#include "Model.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#include <iostream>

Model::Model(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
		: Shape(
			position,
			rotation,
			scale,
			get_geometry(),
			"Model",
			{1.0f, 1.0f, 0.0f},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			}
		)
{

}

Model::~Model() {
	
}


std::vector<GLfloat> Model::get_vertices() {

}

std::vector<GLuint> Model::get_indices() {

}


std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> Model::get_geometry() {
	std::string inputfile = "models/wick.obj";
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

	std::string err;

	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, inputfile.c_str());

	if (!err.empty()) {
		std::cout << err << std::endl;
	}
	if (!ret) {
		return {};
	}

	std::vector<GLfloat> vertices;
	std::vector<GLfloat> normals;
	std::vector<GLuint> indices;

	for (size_t s = 0; s < shapes.size(); s++) {
  		// Loop over faces(polygon)
  		size_t index_offset = 0;
  		for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
    		int fv = shapes[s].mesh.num_face_vertices[f];
			for (size_t v = 0; v < fv; v++) {
				tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
				tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
			  	tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
			  	tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
			  	tinyobj::real_t nx = attrib.normals[3*idx.normal_index+0];
			  	tinyobj::real_t ny = attrib.normals[3*idx.normal_index+1];
			  	tinyobj::real_t nz = attrib.normals[3*idx.normal_index+2];
			  	tinyobj::real_t tx = attrib.texcoords[2*idx.texcoord_index+0];
      			tinyobj::real_t ty = attrib.texcoords[2*idx.texcoord_index+1];
				vertices.push_back(vx);
				vertices.push_back(vy);
				vertices.push_back(vz);
				normals.push_back(nx);
				normals.push_back(ny);
				normals.push_back(nz);
				indices.push_back(3*idx.vertex_index+0);
				indices.push_back(3*idx.vertex_index+1);
				indices.push_back(3*idx.vertex_index+2);
			}
			index_offset += fv;
		}
  	}
	
	return {vertices, normals, indices};
}
