#include "MessageBus.h"

#include "System.h"

MessageBus::MessageBus()
	: _systems{}
{
	
}

void MessageBus::receiveMessage(Message & msg)
{
	switch(msg._type)
	{
		case HELLO:
			{
				auto system = msg._system;
				if (system != nullptr)
				{
					_systems.emplace_back(system);
					Message helloAck {HELLO_ACK};
					system->handleMessage(helloAck);
				}
				else
				{
					cout("Hello from nullptr system");
					exit(1);
				}
			}
			break;

		case DRAW:
			break;

		default:
			break;
	}
	
	// Broadcasting the message to all systems
	// Should optimize it to send it only to the right system
	// (maybe MT this?)
	for (auto system : _systems)
	{
		system->handleMessage(msg);
	}
}

void MessageBus::cout(std::string string) const
{
	std::cout << "0x" << std::hex << std::this_thread::get_id() << " ";
	std::cout << "\033[45m\033[1m";
	std::cout << "[MessageBus]";
	std::cout << "\033[49m\033[0m";
	std::cout << " " << string << std::endl;
}
