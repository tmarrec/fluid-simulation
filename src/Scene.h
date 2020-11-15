#pragma once

#include "glm/ext/scalar_constants.hpp"
#include "glm/gtc/quaternion.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtx/matrix_decompose.hpp"
#include "glm/gtx/string_cast.hpp"
#include "src/SubdivideComponent.h"
#include <GL/gl.h>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define TINYGLTF_USE_CPP14

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#include "tiny_gltf.h"
#pragma GCC diagnostic pop

#include "utils.h"
#include "shapes.h"
#include "ECS.h"
#include "TransformComponent.h"
#include "CameraComponent.h"
#include "DrawableComponent.h"

class Scene
{
public:
	Scene(std::shared_ptr<Renderer> __renderer, std::shared_ptr<ECS_Manager> __ECS_manager, const std::string& sceneFileName)
	: _renderer { __renderer }
	, _ECS_manager { __ECS_manager }
	{
		tinygltf::Model model;
		_readModel(model, sceneFileName);
		_generateEntities(model);
	}


private:
	static void _readModel(tinygltf::Model& __model, const std::string& __fileName)
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
			Shape shape;
			shape = _meshLoop(__model, __model.meshes[__node.mesh]);
			// New entity creation for the mesh
			// TODO SHADER
			auto shader = std::make_shared<Shader>(Shader{"shaders/vert.vert", "shaders/frag.frag"});
			auto& entity = _ECS_manager->addEntity();

			glm::vec3 translation = {0.0f, 0.0f, 0.0f};
			glm::vec3 rotation = {0.0f, 0.0f, 0.0f};
			glm::vec3 scale = {1.0f, 1.0f, 1.0f};
			if (__node.matrix.size() > 0)
			{
				glm::quat qrotation;
				glm::mat4 matrix = glm::make_mat4(&__node.matrix.at(0));
				glm::vec3 skew;
				glm::vec4 perspective;
				glm::decompose(matrix, scale, qrotation, translation, skew, perspective);
				rotation = (glm::eulerAngles(qrotation)/(2*glm::pi<float>()))*glm::vec3{360.0f};
			}
			else
			{
				if (__node.translation.size() > 0)
				{
					translation = glm::make_vec3(&__node.translation.at(0));
				}
				if (__node.rotation.size() > 0)
				{
					rotation = glm::make_vec3(&__node.rotation.at(0));
				}
				if (__node.scale.size() > 0)
				{
					scale = glm::make_vec3(&__node.scale.at(0));
				}
			}
			entity.addComponent<TransformComponent>(translation, rotation, scale);
			entity.addComponent<DrawableComponent>(_renderer, shader, shape.vertices, shape.normals, shape.indices, GL_TRIANGLES);

		}
		for (size_t i = 0; i < __node.children.size(); i++) {
			ASSERT((__node.children[i] >= 0) && (__node.children[i] < __model.nodes.size()), "");
			_nodeLoop(__model, __model.nodes[__node.children[i]]);
		}
	}


	Shape _meshLoop(tinygltf::Model& __model, tinygltf::Mesh& __mesh) const
	{
		std::vector<GLfloat> vertices;
		std::vector<GLfloat> normals;
		std::vector<GLuint> indices;
		for (const auto& primitive : __mesh.primitives)
		{
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
					vertices.insert(vertices.begin(), verticesBuffer, verticesBuffer+accessor.count*3);
				}
				else if (attribute.first.compare("NORMAL") == 0)
				{
					GLfloat* positionBuffer = reinterpret_cast<GLfloat*>(buffer);
					normals.insert(normals.begin(), positionBuffer, positionBuffer+accessor.count*3);
				}
			}
			
			// Indices
			if (primitive.indices >= 0)
			{
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
			}
			else
			{
				ASSERT(vertices.size()%3==0, "vertices must be power of 3");
				indices.resize(vertices.size()/3);
				std::iota(indices.begin(), indices.end(), 0);
			}
			
 		} 

		return {vertices, normals, indices};
	}

	template<typename T>
	void _getIndicesFromBuffer(std::vector<GLuint>& __indices, unsigned char* __buffer, const std::uint64_t __count) const
	{
		T* indicesBuffer = reinterpret_cast<T*>(__buffer);
		__indices.reserve(__count);
		for (std::uint64_t i = 0; i < __count; ++i)
		{
			__indices.emplace_back(static_cast<GLuint>(indicesBuffer[i]));
		}
	}
	
	std::shared_ptr<Renderer> _renderer;
	std::shared_ptr<ECS_Manager> _ECS_manager;
};
