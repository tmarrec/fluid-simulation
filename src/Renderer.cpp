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

void Renderer::initGl(int width, int height) const
{
	glEnable(GL_DEPTH_TEST);
	glViewport(0, 0, width, height);

	cout(std::string("Renderer       : ")+reinterpret_cast<char const*>(glGetString(GL_RENDERER)));
	cout(std::string("Vendor         : ")+reinterpret_cast<char const*>(glGetString(GL_VENDOR)));
	cout(std::string("Version        : ")+reinterpret_cast<char const*>(glGetString(GL_VERSION)));
	cout(std::string("GLSL Version   : ")+reinterpret_cast<char const*>(glGetString(GL_SHADING_LANGUAGE_VERSION)));
}

void Renderer::resizeGl(int width, int height) const
{
	glViewport(0, 0, width, height);
}

void Renderer::initDrawable(std::uint64_t componentID, Shape shape)
{
	glGenBuffers(1, &_GLObjects[componentID].VBO);
	glGenBuffers(1, &_GLObjects[componentID].NBO);
	glGenBuffers(1, &_GLObjects[componentID].EBO);

	glGenVertexArrays(1, &_GLObjects[componentID].VAO);

	glBindVertexArray(_GLObjects[componentID].VAO);

		glBindBuffer(GL_ARRAY_BUFFER, _GLObjects[componentID].VBO);
		glBufferData(
			GL_ARRAY_BUFFER,
			shape.vertices.size()*sizeof(GLfloat),
			shape.vertices.data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, _GLObjects[componentID].NBO);
		glBufferData(
			GL_ARRAY_BUFFER,
			shape.normals.size()*sizeof(GLfloat),
			shape.normals.data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _GLObjects[componentID].EBO);
		glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			shape.indices.size()*sizeof(GLuint),
			shape.indices.data(),
			GL_STATIC_DRAW
		);

	glBindVertexArray(0);
}

void Renderer::freeDrawable(std::uint64_t componentID, Shape shape)
{
	glDeleteVertexArrays(1, &_GLObjects[componentID].VAO);	
	glDeleteBuffers(1, &_GLObjects[componentID].VBO);
	glDeleteBuffers(1, &_GLObjects[componentID].NBO);
	glDeleteBuffers(1, &_GLObjects[componentID].EBO);
}

void Renderer::draw(std::uint64_t componentID, Shape shape)
{
	cout("DRAW");
	cout(std::to_string(componentID));
	glClearColor(1.0f, 1.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glBindVertexArray(_GLObjects[componentID].VAO);
	glDrawElements(GL_TRIANGLES, shape.indices.size(), GL_UNSIGNED_INT, nullptr);
	glBindVertexArray(0);
}

void Renderer::handleMessage(Message & msg)
{
	switch(msg._type)
	{
		case HELLO_ACK:
			cout("Loaded in the \033[45m\033[1m[MessageBus]\033[49m\033[0m");
			break;

		case INIT_GL:
			initGl(msg._width, msg._height);
			break;

		case RESIZE_GL:
			resizeGl(msg._width, msg._height);
			break;

		case INIT_DRAWABLE:
			initDrawable(msg._componentID, msg._shape);
			break;

		case FREE_DRAWABLE:
			freeDrawable(msg._componentID, msg._shape);
			break;

		case DRAW:
			draw(msg._componentID, msg._shape);
			break;

		default:
			break;
	}
}

