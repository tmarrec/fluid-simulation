#pragma once

#include <vector>
#include <string>

#include "Message.h"

class System;

class MessageBus
{
public:
	MessageBus();
	void receiveMessage(Message & msg);
	void cout(std::string string) const;

private:
	std::vector<System *> _systems;

};

