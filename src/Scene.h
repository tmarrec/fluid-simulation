#pragma once

#include <GL/gl.h>
#include <cstdint>
#include <iostream>

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define TINYGLTF_USE_CPP14

#include "tinygltf/tiny_gltf.h"

#include "utils.h"

class Scene
{
public:
	Scene()
	{
		tinygltf::Model model;
		tinygltf::TinyGLTF loader;
		std::string err;
		std::string warn;
		bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, "scene.gltf");
		if (!warn.empty()) {
			printf("Warn: %s\n", warn.c_str());
		}
		if (!err.empty()) {
			printf("Err: %s\n", err.c_str());
		}
		if (!ret) {
			printf("Failed to parse glTF\n");
			exit(-1);
		}

		std::map<int, GLuint> vbos;

		const tinygltf::Scene &scene = model.scenes[model.defaultScene];
		for (size_t i = 0; i < scene.nodes.size(); ++i) {
			assert((scene.nodes[i] >= 0) && (scene.nodes[i] < model.nodes.size()));
			bindModelNodes(vbos, model, model.nodes[scene.nodes[i]]);
		}

	}

	void bindModelNodes(std::map<int, GLuint> vbos, tinygltf::Model &model, tinygltf::Node &node)
	{
		if ((node.mesh >= 0) && (node.mesh < model.meshes.size())) {
			bindMesh(vbos, model, model.meshes[node.mesh]);
		}
		for (size_t i = 0; i < node.children.size(); i++) {
			assert((node.children[i] >= 0) && (node.children[i] < model.nodes.size()));
			bindModelNodes(vbos, model, model.nodes[node.children[i]]);
		}
	}

	std::map<int, GLuint> bindMesh(std::map<int, GLuint> vbos, tinygltf::Model &model, tinygltf::Mesh &mesh)
	{
		for (const auto& primitive : mesh.primitives)
		{
			std::cout << "new entity" << std::endl;
			std::vector<float> vertices;
    		for (const auto& attribute : primitive.attributes)
			{
        		const auto& accessor = model.accessors[attribute.second];
        		if (attribute.first.compare("POSITION") == 0)
				{
					// Get vertex buffer
					const auto& bufferView = model.bufferViews[accessor.bufferView];
            		auto& vertexArray = model.buffers[bufferView.buffer];
					unsigned char* buffer = &vertexArray.data.at(0) + bufferView.byteOffset + accessor.byteOffset;

					// Cast the vertex buffer to clean float vector
					float* verticesBuffer = reinterpret_cast<float*>(buffer);
					vertices.insert(vertices.cbegin(), verticesBuffer, verticesBuffer+accessor.count);
        		}
    		}
			for (auto v : vertices)
			{
				std::cout << v << " ";
			}
			std::cout << std::endl;
 		} 
		return vbos;
	}

private:
};
