#pragma once
#include <vulkan/vulkan.h>
#include <utility>
#include <cstdint>
#include <vector>
#include <optional>

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

struct QueueFamilyIndices
{
	std::optional<std::uint32_t> graphics;

	bool isComplete()
	{
		return graphics.has_value();
	}
};

class Renderer
{
public:
	void init(std::shared_ptr<Window> window);
	~Renderer();

private:
	void createInstance();
	bool checkValidationLayerSupport();
	std::vector<const char*> getVkRequiredExtensions();
	void setupDebugMessenger();
	void populateDebugMessengerCreateInfo(VkDebugUtilsMessengerCreateInfoEXT& createInfo);
	void destroyDebugMessenger();
	void pickPhysicalDevice();
	bool isPhysicalDeviceSuitable(VkPhysicalDevice device);
	QueueFamilyIndices findQueueFamilies(VkPhysicalDevice device);

	VkInstance _vkInstance = nullptr;
	std::shared_ptr<Window> _window = nullptr;
	VkDebugUtilsMessengerEXT _debugMessenger = nullptr;
	VkPhysicalDevice _physicalDevice = VK_NULL_HANDLE;
};


static VKAPI_ATTR VkBool32 VKAPI_CALL vkDebugCallback
(
	VkDebugUtilsMessageSeverityFlagBitsEXT messageSeverity,
	VkDebugUtilsMessageTypeFlagsEXT messageType,
	const VkDebugUtilsMessengerCallbackDataEXT* pCallbackData,
	void* pUserData
)
{
	switch (messageSeverity)
	{
		case VK_DEBUG_UTILS_MESSAGE_SEVERITY_INFO_BIT_EXT:
			INFO(pCallbackData->pMessage);
			break;
		case VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT:
			WARNING(pCallbackData->pMessage);
			break;
		case VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT:
			ERROR(pCallbackData->pMessage);
			break;
		default:
			break;
	}
	return VK_FALSE;
}
