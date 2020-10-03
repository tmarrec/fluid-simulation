#pragma once

#include "Message.h"

#include <iostream>
#include <thread>

class MessageBus;

class System
{
public:
	System(MessageBus & messageBus);
	void handleMessage(Message & msg) const;
	void postMessage(Message & msg) const;
	MessageBus & __messageBus(); // Cannot name it messageBus()
	virtual void cout(std::string string) const = 0;

private:
	MessageBus & _messageBus;
};

