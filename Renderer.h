#pragma once
#include <vulkan/vulkan.h>
#include <utility>
#include <cstdint>
#include <vector>

#include "Window.h"

class Window;

class Renderer
{
public:
	void init(std::shared_ptr<Window> window);
	~Renderer();

private:
	void createInstance();

	VkInstance _vkInstance = nullptr;
	std::shared_ptr<Window> _window = nullptr;
};
