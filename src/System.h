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

using MsgBus_ptr = std::shared_ptr<MessageBus>;

class System
{
public:
	System(MsgBus_ptr messageBus);
	virtual void handleMessage(Message & msg) = 0;
	void postMessage(Message & msg) const;
	MsgBus_ptr __messageBus(); // Cannot name it messageBus()
	virtual void cout(std::string string) const = 0;

private:
	MsgBus_ptr _messageBus;
};

