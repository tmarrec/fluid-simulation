#include "Renderer.h"

Renderer::Renderer(MessageBus &	messageBus)
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

void Renderer::draw() const
{
	cout("DRAW");
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
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

void Renderer::handleMessage(Message & msg) const
{
	switch(msg.type())
	{
		case HELLO_ACK:
			cout("Loaded in the \033[45m\033[1m[MessageBus]\033[49m\033[0m");
			break;

		case INIT_GL:
			initGl(msg.width(), msg.height());
			break;

		case RESIZE_GL:
			resizeGl(msg.width(), msg.height());
			break;

		case DRAW:
			draw();
			break;

		default:
			break;
	}
}

