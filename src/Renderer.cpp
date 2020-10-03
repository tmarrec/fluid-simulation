#include "Renderer.h"

Renderer::Renderer(MessageBus &	messageBus)
	: System{messageBus}
{

}

void Renderer::exec()
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
