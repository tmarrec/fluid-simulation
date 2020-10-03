#pragma once

#include "System.h"

class Renderer : public System
{
public:
	Renderer(MessageBus & messageBus);
	void exec();
	void cout(std::string string) const override;

private:

};
