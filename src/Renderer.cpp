#include "Renderer.h"


Renderer::Renderer(MsgBus_ptr messageBus)
	: System{messageBus}
{
	Message helloMsg (HELLO, this);
	postMessage(helloMsg);
}

void Renderer::cout(std::string string) const
{
	std::cout << "0x" << std::hex << std::this_thread::get_id() << " ";
	std::cout << "  \033[44m\033[1m";
	std::cout << "[Renderer]";
	std::cout << "\033[49m\033[0m";
	std::cout << " " << string << std::endl;
}

void Renderer::initGl(Message & msg)
{
	_glContext = msg._glContext;
	_glWidget = msg._glWidget;

	//_glWidget->makeCurrent();
	glEnable(GL_DEPTH_TEST);
	glViewport(0, 0, msg._width, msg._height);
	/*
	cout(std::string("Renderer       : ")+reinterpret_cast<char const*>(glGetString(GL_RENDERER)));
	cout(std::string("Vendor         : ")+reinterpret_cast<char const*>(glGetString(GL_VENDOR)));
	cout(std::string("Version        : ")+reinterpret_cast<char const*>(glGetString(GL_VERSION)));
	cout(std::string("GLSL Version   : ")+reinterpret_cast<char const*>(glGetString(GL_SHADING_LANGUAGE_VERSION)));
	*/
	//_glWidget->doneCurrent();
	cout("OpenGL initialized");
}

void Renderer::resizeGl(int width, int height) const
{
	glViewport(0, 0, width, height);
}

void Renderer::initDrawable(Message & msg)
{
	// not sure about this
	_glWidget->makeCurrent();
	//_glContext->doneCurrent();

	//std::cout << "isvalid: " << _glContext->isValid() << std::endl;
	std::cout << "start init" << std::endl;
	glGenBuffers(1, &msg._VBO);
	glGenBuffers(1, &msg._NBO);
	glGenBuffers(1, &msg._EBO);

	glGenVertexArrays(1, &msg._VAO);

	std::cout << msg._vertices->size() << std::endl;
	std::cout << msg._vertices->data() << std::endl;

	glBindVertexArray(msg._VAO);
		glBindBuffer(GL_ARRAY_BUFFER, msg._VBO);
		glBufferData(
			GL_ARRAY_BUFFER,
			msg._vertices->size()*sizeof(GLfloat),
			msg._vertices->data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, msg._NBO);
		glBufferData(
			GL_ARRAY_BUFFER,
			msg._normals->size()*sizeof(GLfloat),
			msg._normals->data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, msg._EBO);
		glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			msg._indices->size()*sizeof(GLuint),
			msg._indices->data(),
			GL_STATIC_DRAW
		);

	glBindVertexArray(0);
	std::cout << "end init" << std::endl;

	_glWidget->doneCurrent();
}

void Renderer::freeDrawable(Message & msg)
{
	/* THIS AS NO SENSE IF IT'S NOT A POINTER */
	/*
	glDeleteVertexArrays(1, &_componentAttributes[componentID]._GLObjects.VAO);	
	glDeleteBuffers(1, &_componentAttributes[componentID]._GLObjects.VBO);
	glDeleteBuffers(1, &_componentAttributes[componentID]._GLObjects.NBO);
	glDeleteBuffers(1, &_componentAttributes[componentID]._GLObjects.EBO);
	*/
}

glm::mat4 Renderer::_getModel(Message & msg)
{
	glm::mat4 model {1.0f};
	model = glm::translate(model, msg._position);
	model = glm::rotate(model, glm::radians(msg._rotation.x), glm::vec3(1.0f, 0.0f, 0.0f));
	model = glm::rotate(model, glm::radians(msg._rotation.y), glm::vec3(0.0f, 1.0f, 0.0f));
	model = glm::rotate(model, glm::radians(msg._rotation.z), glm::vec3(0.0f, 0.0f, 1.0f));
	model = glm::scale(model, glm::vec3{msg._scale});

	return model;
}

void Renderer::_useShader(Message & msg)
{
	auto shader = msg._shader;
	auto color = msg._color;
	shader->use();
	shader->set_3f("_object_color", color);

	shader->set_mat4("model", _getModel(msg));
	shader->set_mat4("view", msg._view);
	shader->set_mat4("projection", msg._projection);
	uint64_t time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	float time_sin = sin(time/50000);
	shader->set_1f("time", time_sin);

	shader->set_3f("_view_pos", msg._view[3]);

	// TODO add ligths 	
	shader->set_1i("_light_nb", 0);

	/*
	shader->set_1i("_light_nb", lights.size());
	// Envoie les uniforms pour toutes les lumieres
	for (size_t i = 0; i < lights.size(); ++i) {
		auto light = std::static_pointer_cast<Light>(lights[i]);
		auto temp = std::string("_point_lights[") + std::to_string(i) + "].position";
		shader->set_3f(temp.c_str(), light->position());
		temp = std::string("_point_lights[") + std::to_string(i) + "].color";
		shader->set_3f(temp.c_str(), light->color());
		temp = std::string("_point_lights[") + std::to_string(i) + "].intensity";
		shader->set_1f(temp.c_str(), light->intensity());
		temp = std::string("_point_lights[") + std::to_string(i) + "].constant";
		shader->set_1f(temp.c_str(), 1.0f);
		temp = std::string("_point_lights[") + std::to_string(i) + "].linear";
		shader->set_1f(temp.c_str(), 0.0014f);
		temp = std::string("_point_lights[") + std::to_string(i) + "].quadratic";
		shader->set_1f(temp.c_str(), 0.000007f);
	}
	*/
}

void Renderer::draw(Message & msg)
{
	//_glContext->doneCurrent();
	//_glWidget->makeCurrent();

	glClearColor(0.8f, 0.8f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	_useShader(msg);
	
	std::cout << msg._VAO << std::endl;
	std::cout << msg._indices->size() << std::endl;

	glBindVertexArray(msg._VAO);
	std::cout << "1" << std::endl;
	glDrawElements(GL_TRIANGLES, msg._indices->size(), GL_UNSIGNED_INT, nullptr);
	std::cout << "2" << std::endl;
	glBindVertexArray(0);

	glFinish();
	//_glWidget->doneCurrent();
}

void Renderer::handleMessage(Message & msg)
{
	switch(msg._type)
	{
		case HELLO_ACK:
			cout("Loaded in the \033[45m\033[1m[MessageBus]\033[49m\033[0m");
			break;

		case INIT_GL:
			initGl(msg);
			break;

		case RESIZE_GL:
			resizeGl(msg._width, msg._height);
			break;

		case INIT_DRAWABLE:
			initDrawable(msg);
			break;

		case FREE_DRAWABLE:
			freeDrawable(msg);
			break;

		case DRAW:
			cout("Draw");
			draw(msg);
			break;

		default:
			break;
	}
}

