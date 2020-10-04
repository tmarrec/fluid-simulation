#include "System.h" 

#include "MessageBus.h"

System::System(MsgBus_ptr messageBus)
: _messageBus{messageBus}
{

}

void System::postMessage(Message & msg) const
{
	_messageBus->receiveMessage(msg);
}

MsgBus_ptr System::__messageBus()
{
	return _messageBus;
}
