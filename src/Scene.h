#pragma once

#include <GL/gl.h>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define TINYGLTF_USE_CPP14

#include "tinygltf/tiny_gltf.h"

#include "utils.h"
#include "ECS.h"
#include "TransformComponent.h"
#include "CameraComponent.h"
#include "DrawableComponent.h"

class Scene
{
public:
	Scene(std::shared_ptr<Renderer> __renderer, std::shared_ptr<ECS_Manager> __ECS_manager, const std::string sceneFileName)
	: _renderer { __renderer }
	, _ECS_manager { __ECS_manager }
	{
		tinygltf::Model model;
		_readModel(model, sceneFileName);
		_generateEntities(model);
	}


private:
	void _readModel(tinygltf::Model& __model, const std::string __fileName) const
	{
		tinygltf::TinyGLTF loader;
		std::string err;
		std::string warn;
		bool ret = loader.LoadASCIIFromFile(&__model, &err, &warn, __fileName);
		if (!warn.empty())
		{
			WARNING(warn.c_str());
		}
		if (!err.empty()) 
		{
			ERROR(err.c_str());
		}
		if (!ret)
		{
			ERROR("Failed to parse glTF file");
		}
	}

	void _generateEntities(tinygltf::Model& __model)
	{
		const tinygltf::Scene &scene = __model.scenes[__model.defaultScene];
		// TODO handle multiples scenes
		for (size_t i = 0; i < scene.nodes.size(); ++i)
		{
			ASSERT((scene.nodes[i] >= 0) && (scene.nodes[i] < __model.nodes.size()), "");
			_nodeLoop(__model, __model.nodes[scene.nodes[i]]);
		}
	}

	void _nodeLoop(tinygltf::Model& __model, tinygltf::Node& __node)
	{
		if ((__node.mesh >= 0) && (__node.mesh < (int)__model.meshes.size())) {
			_meshLoop(__model, __model.meshes[__node.mesh]);
		}
		for (size_t i = 0; i < __node.children.size(); i++) {
			ASSERT((node.children[i] >= 0) && (node.children[i] < model.nodes.size()), "");
			_nodeLoop(__model, __model.nodes[__node.children[i]]);
		}
	}


	void _meshLoop(tinygltf::Model& __model, tinygltf::Mesh& __mesh)
	{
		for (const auto& primitive : __mesh.primitives)
		{
			std::vector<GLfloat> vertices;
			std::vector<GLfloat> normals;
			std::vector<GLuint> indices;
    		for (const auto& attribute : primitive.attributes)
			{
        		const auto& accessor = __model.accessors[attribute.second];
				const auto& bufferView = __model.bufferViews[accessor.bufferView];
            	auto& bufferArray = __model.buffers[bufferView.buffer];
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
			const auto& indicesAccessor = __model.accessors[primitive.indices];
			const auto& type = indicesAccessor.componentType;
			auto& indicesArray = __model.buffers[__model.bufferViews[indicesAccessor.bufferView].buffer];
			const auto& indicesBufferView = __model.bufferViews[indicesAccessor.bufferView];
			unsigned char* buffer = &indicesArray.data.at(0) + indicesBufferView.byteOffset + indicesAccessor.byteOffset;
			switch (type)
			{
				case TINYGLTF_COMPONENT_TYPE_BYTE:
					_getIndicesFromBuffer<std::int8_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:
					_getIndicesFromBuffer<std::uint8_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_SHORT:
					_getIndicesFromBuffer<std::int16_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT:
					_getIndicesFromBuffer<std::uint16_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_INT:
					_getIndicesFromBuffer<std::int32_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT:
					_getIndicesFromBuffer<std::uint32_t>(indices, buffer, indicesAccessor.count);
					break;
				case TINYGLTF_COMPONENT_TYPE_FLOAT:
				case TINYGLTF_COMPONENT_TYPE_DOUBLE:
					ERROR("Indices componentType can't be a float or a double");
					break;
			}
			
			std::cout << "new entity" << std::endl;
			auto shader = std::make_shared<Shader>(Shader{"shaders/vert.vert", "shaders/frag.frag"});
			auto& entity = _ECS_manager->addEntity();
			entity.addComponent<TransformComponent>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{20.0f, 20.0f, 20.0f});
			entity.addComponent<DrawableComponent>(_renderer, shader, vertices, normals, indices, GL_TRIANGLES);
 		} 
	}

	template<typename T>
	void _getIndicesFromBuffer(std::vector<GLuint>& __indices, unsigned char* __buffer, const std::uint64_t __count) const
	{
		T* indicesBuffer = reinterpret_cast<T*>(__buffer);
		for (std::uint64_t i = 0; i < __count; ++i)
		{
			__indices.emplace_back(static_cast<GLuint>(indicesBuffer[i]));
		}
	}
	
	std::shared_ptr<Renderer> _renderer;
	std::shared_ptr<ECS_Manager> _ECS_manager;
};
