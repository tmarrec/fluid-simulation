#include "System.h" 

#include "MessageBus.h"

System::System(MessageBus & messageBus)
: _messageBus{messageBus}
{

}

void System::postMessage(Message & msg) const
{
	_messageBus.receiveMessage(msg);
}

MessageBus & System::__messageBus()
{
	return _messageBus;
}
