#include "Renderer.h"

void Renderer::init(std::shared_ptr<Window> window)
{
	_window = window;
	createInstance();
}

Renderer::~Renderer()
{
	vkDestroyInstance(_vkInstance, nullptr);
}

void Renderer::createInstance()
{
	ASSERT(_window, "_window pointer should not be nullptr");

	if (enableValidationLayers && !checkValidationLayerSupport())
	{
		ERROR("Validation layers requested but not available.");
	}

	VkApplicationInfo appInfo{};
	appInfo.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
	appInfo.pNext = nullptr;
	appInfo.pApplicationName = "cowboy-engine";
	appInfo.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
	appInfo.pEngineName = "cowboy-engine";
	appInfo.engineVersion = VK_MAKE_VERSION(0, 1, 0);
	appInfo.apiVersion = VK_API_VERSION_1_2;

	VkInstanceCreateInfo createInfo{};
	createInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
	createInfo.pNext = nullptr;
	createInfo.flags = 0;
	createInfo.pApplicationInfo = &appInfo;
	if (enableValidationLayers)
	{
		createInfo.enabledLayerCount = static_cast<std::uint32_t>(vkValidationLayers.size());
		createInfo.ppEnabledLayerNames = vkValidationLayers.data();
	}
	else
	{
		createInfo.enabledLayerCount = 0;
		createInfo.ppEnabledLayerNames = nullptr;
	}
	std::pair<const char**, std::uint32_t> windowRequiredInstanceExtensions = _window->windowGetRequiredInstanceExtensions();
	createInfo.enabledExtensionCount = windowRequiredInstanceExtensions.second;
	createInfo.ppEnabledExtensionNames = windowRequiredInstanceExtensions.first;

	if (vkCreateInstance(&createInfo, nullptr, &_vkInstance))
	{
		ERROR("Cannot create Vulkan instance.");
	}
}

bool Renderer::checkValidationLayerSupport()
{
	std::uint32_t vkLayerCount;
	if (vkEnumerateInstanceLayerProperties(&vkLayerCount, nullptr))
	{
		ERROR("Cannot enumerate instance layer properties.");
	}

	std::vector<VkLayerProperties> vkAvailableLayers{ vkLayerCount };
	if (vkEnumerateInstanceLayerProperties(&vkLayerCount, vkAvailableLayers.data()))
	{
		ERROR("Cannot enumerate instance layer properties.");
	}

	// Search if all of the layers in vkValidationLayers exist in vkAvailableLayers
	for (const char* layerName : vkValidationLayers)
	{
		bool layerFound = false;
		for (const VkLayerProperties& layerProperties : vkAvailableLayers)
		{
			if (strcmp(layerName, layerProperties.layerName))
			{
				layerFound = true;
				break;
			}
		}
		if (!layerFound)
		{
			return false;
		}
	}
	return true;
}
