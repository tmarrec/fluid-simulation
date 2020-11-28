#pragma once
#include <vulkan/vulkan.h>
#include <utility>
#include <cstdint>
#include <vector>

#include "Window.h"

const std::vector<const char*> vkValidationLayers =
{
	"VK_LAYER_KHRONOS_validation"
};

#ifdef NDEBUG
const bool enableValidationLayers = false;
#else
const bool enableValidationLayers = true;
#endif

class Renderer
{
public:
	void init(std::shared_ptr<Window> window);
	~Renderer();

private:
	void createInstance();
	bool checkValidationLayerSupport();

	VkInstance _vkInstance = nullptr;
	std::shared_ptr<Window> _window = nullptr;
};
