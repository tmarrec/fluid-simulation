#pragma once

#include "Window.h"

class Renderer
{
public:
	void init(std::shared_ptr<Window> window);
    void pass();
	~Renderer();

private:

	std::shared_ptr<Window> _window = nullptr;
};

