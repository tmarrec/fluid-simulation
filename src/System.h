#pragma once

#include "Message.h"

#ifdef __APPLE__
	#include <OpenGL/gl.h>
#else
	#include <GL/gl.h>
#endif

#include <iostream>
#include <thread>

class MessageBus;

class System
{
public:
	System(MessageBus & messageBus);
	virtual void handleMessage(Message & msg) const = 0;
	void postMessage(Message & msg) const;
	MessageBus & __messageBus(); // Cannot name it messageBus()
	virtual void cout(std::string string) const = 0;

private:
	MessageBus & _messageBus;
};

