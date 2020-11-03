#pragma once

#include <GL/gl.h>
#include <cstddef>
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

	template<typename T>
	void getIndicesFromBuffer(std::vector<GLuint>& indices, unsigned char* buffer, const size_t count)
	{
		T* indicesBuffer = reinterpret_cast<T*>(buffer);
		for (std::uint64_t i = 0; i < count; ++i)
		{
			indices.emplace_back(static_cast<GLuint>(indicesBuffer[i]));
		}
	}

	std::map<int, GLuint> bindMesh(std::map<int, GLuint> vbos, tinygltf::Model &model, tinygltf::Mesh &mesh)
	{
		for (const auto& primitive : mesh.primitives)
		{
			std::cout << "new entity" << std::endl;
			std::vector<GLfloat> vertices;
			std::vector<GLfloat> normals;
			std::vector<GLuint> indices;
    		for (const auto& attribute : primitive.attributes)
			{
        		const auto& accessor = model.accessors[attribute.second];
				const auto& bufferView = model.bufferViews[accessor.bufferView];
            	auto& bufferArray = model.buffers[bufferView.buffer];
				unsigned char* buffer = &bufferArray.data.at(0)+bufferView.byteOffset+accessor.byteOffset;
        		if (attribute.first.compare("POSITION") == 0)
				{
					// Cast the vertex array to a clean float vector
					GLfloat* verticesBuffer = reinterpret_cast<GLfloat*>(buffer);
					vertices.insert(vertices.begin(), verticesBuffer, verticesBuffer+accessor.count);
        		}
				else if (attribute.first.compare("NORMAL") == 0)
				{
					GLfloat* positionBuffer = reinterpret_cast<GLfloat*>(buffer);
					normals.insert(normals.begin(), positionBuffer, positionBuffer+accessor.count);
				}
    		}
			
			// Indices
			const auto& indicesAccessor = model.accessors[primitive.indices];
			auto& type = indicesAccessor.componentType;
			std::cout << type << std::endl;
			auto& indicesArray = model.buffers[model.bufferViews[indicesAccessor.bufferView].buffer];
			auto& indicesBufferView = model.bufferViews[indicesAccessor.bufferView];
			unsigned char* buffer = &indicesArray.data.at(0) + indicesBufferView.byteOffset + indicesAccessor.byteOffset;
			switch (type)
			{
				case TINYGLTF_COMPONENT_TYPE_BYTE:
					getIndicesFromBuffer<std::int8_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:
					getIndicesFromBuffer<std::uint8_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_SHORT:
					getIndicesFromBuffer<std::int16_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT:
					getIndicesFromBuffer<std::uint16_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_INT:
					getIndicesFromBuffer<std::int32_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT:
					getIndicesFromBuffer<std::uint32_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_FLOAT:
				case TINYGLTF_COMPONENT_TYPE_DOUBLE:
					ERROR("Indices componentType can't be a float or a double");
					break;
			}

			std::cout << "Vertices : " << std::endl;
			for (auto v : vertices)
			{
				std::cout << v << " ";
			}
			std::cout << std::endl;
			std::cout << "Normals : " << std::endl;
			for (auto v : normals)
			{
				std::cout << v << " ";
			}
			std::cout << std::endl;
			std::cout << "Indices : " << std::endl;
			for (auto v : indices)
			{
				std::cout << v << " ";
			}
			std::cout << std::endl;
 		} 
		return vbos;
	}

private:
};
