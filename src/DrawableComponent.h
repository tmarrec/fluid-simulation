#pragma once

#include "ECS.h"
#include "Shader.h"

class DrawableComponent : public Component, public System
{
private:
	std::vector<GLfloat> _vertices;
	std::vector<GLfloat> _normals;
	std::vector<GLuint>  _indices;

	GLuint _VAO;
	GLuint _VBO;
	GLuint _NBO;
	GLuint _EBO;

	glm::vec3 _color = {1.0f, 0.5f, 0.5f};
	std::unique_ptr<Shader> _shader;


public:
	DrawableComponent(std::string vertPath, std::string fragPath,
		std::vector<GLfloat> vertices, std::vector<GLfloat> normals,
		std::vector<GLuint> indices, MsgBus_ptr messageBus)
		: System{messageBus}
	{
		_vertices = vertices;
		_normals = normals;
		_indices = indices;


		auto shader = new Shader{vertPath.c_str(), fragPath.c_str()};
		std::unique_ptr<Shader> uPtr {shader};
		_shader = std::move(uPtr);

		Message helloMsg {HELLO, this};
		postMessage(helloMsg);
	}

	void cout(std::string string) const
	{
		std::cout << "0x" << std::hex << std::this_thread::get_id() << " ";
		std::cout << "  \033[100m\033[1m";
		std::cout << "[Entity " << getComponentID() << "]";
		std::cout << "\033[49m\033[0m";
		std::cout << " " << string << std::endl;
	}
		
	void handleMessage(Message & msg)
	{
		switch(msg._type)
		{
			case HELLO_ACK:
				cout("Loaded in the \033[45m\033[1m[MessageBus]\033[49m\033[0m");
				break;

			case ASK_ENTITIES_DRAW:
				draw();
				break;
	
			default:
				break;
		}
	}

	void init() override
	{
		Shape shape;
		shape.vertices = _vertices;
		shape.normals = _normals;
		shape.indices = _indices;
		Message initMsg {INIT_DRAWABLE, getComponentID(), shape};
		postMessage(initMsg);
	}

	void draw() override
	{
		Shape shape;
		shape.vertices = _vertices;
		shape.normals = _normals;
		shape.indices = _indices;
		Message drawMsg {DRAW, getComponentID(), shape};
		postMessage(drawMsg);

		//TODO Not sure about this..
		
		/*
		_shader->use();
		_shader->set_3f("_object_color", _color);

		_shader->set_mat4("model", get_model());
		_shader->set_mat4("view", view);
		_shader->set_mat4("projection", projection);

		uint64_t time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		float time_sin = sin(time/50000);
		_shader->set_1f("time", time_sin);

		_shader->set_3f("_view_pos", view[3]);
		_shader->set_1i("_light_nb", lights.size());

		// Envoie les uniforms pour toutes les lumieres
		for (size_t i = 0; i < lights.size(); ++i) {
			auto light = std::static_pointer_cast<Light>(lights[i]);
			auto temp = std::string("_point_lights[") + std::to_string(i) + "].position";
			_shader->set_3f(temp.c_str(), light->position());
			temp = std::string("_point_lights[") + std::to_string(i) + "].color";
			_shader->set_3f(temp.c_str(), light->color());
			temp = std::string("_point_lights[") + std::to_string(i) + "].intensity";
			_shader->set_1f(temp.c_str(), light->intensity());
			temp = std::string("_point_lights[") + std::to_string(i) + "].constant";
			_shader->set_1f(temp.c_str(), 1.0f);
			temp = std::string("_point_lights[") + std::to_string(i) + "].linear";
			_shader->set_1f(temp.c_str(), 0.0014f);
			temp = std::string("_point_lights[") + std::to_string(i) + "].quadratic";
			_shader->set_1f(temp.c_str(), 0.000007f);
		}

		//
		//
		*/
	}

	void update() override
	{
	}

	~DrawableComponent() override
	{
		Shape shape;
		shape.vertices = _vertices;
		shape.normals = _normals;
		shape.indices = _indices;
		Message freeMsg {FREE_DRAWABLE, getComponentID(), shape};
		postMessage(freeMsg);
	}
};
