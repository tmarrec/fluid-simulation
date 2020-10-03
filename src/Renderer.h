#pragma once

#include "System.h"

class Renderer : public System
{
public:
	Renderer(MessageBus & messageBus);
	void cout(std::string string) const override;
	void handleMessage(Message & msg) const override;

	void draw() const;
	void initGl(int width, int height) const;
	void resizeGl(int width, int height) const;

private:

};
