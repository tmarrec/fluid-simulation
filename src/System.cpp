#include "System.h" 

#include "MessageBus.h"

System::System(MessageBus & messageBus)
: _messageBus{messageBus}
{

}

void System::handleMessage(Message & msg) const
{
	switch(msg.type())
	{
		case HELLO_ACK:
			cout("Loaded in the \033[45m\033[1m[MessageBus]\033[49m\033[0m");
			break;			
	}
}

void System::postMessage(Message & msg) const
{
	_messageBus.receiveMessage(msg);
}

MessageBus & System::__messageBus()
{
	return _messageBus;
}
