#pragma once
#include <utility>
#include <cstdint>
#include <vector>
#include <optional>
#include <set>

#include "Window.h"

class Renderer
{
public:
	void init(std::shared_ptr<Window> window);
	~Renderer();

private:

	std::shared_ptr<Window> _window = nullptr;
};

