#include "Renderer.h"

void Renderer::init(std::shared_ptr<Window> window)
{
	_window = window;
	createInstance();
	setupDebugMessenger();
	createSurface();
	pickPhysicalDevice();
	createLogicalDevice();
}

Renderer::~Renderer()
{
	destroyDebugMessenger();
	vkDestroyDevice(_vkDevice, nullptr);
	vkDestroySurfaceKHR(_vkInstance, _vkSurface, nullptr);
	vkDestroyInstance(_vkInstance, nullptr);
}

void Renderer::pickPhysicalDevice()
{
	std::uint32_t deviceCount = 0;
	vkEnumeratePhysicalDevices(_vkInstance, &deviceCount, nullptr);
	if (deviceCount == 0)
	{
		ERROR("Failed to find GPUs with Vulkan support.");
	}
	std::vector<VkPhysicalDevice> devices{ deviceCount };
	vkEnumeratePhysicalDevices(_vkInstance, &deviceCount, devices.data());
	
	for (const auto& device : devices)
	{
		if (isPhysicalDeviceSuitable(device))
		{
			_vkPhysicalDevice = device;
			break;
		}
	}
	
	if (_vkPhysicalDevice == VK_NULL_HANDLE)
	{
		ERROR("Failed to find suitable GPU.");
	}
}

bool Renderer::isPhysicalDeviceSuitable(VkPhysicalDevice device)
{
	VkPhysicalDeviceProperties deviceProperties;
	vkGetPhysicalDeviceProperties(device, &deviceProperties);

	VkPhysicalDeviceFeatures deviceFeatures;
	vkGetPhysicalDeviceFeatures(device, &deviceFeatures);

	QueueFamilyIndices queueFamilyIndices = findQueueFamilies(device);

	// Temporary using only discrete GPU
	if (deviceProperties.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU && deviceFeatures.geometryShader && queueFamilyIndices.isComplete())
	{
		INFO("Found suitable GPU: " << deviceProperties.deviceName);
		return true;
	}
	return false;
}

QueueFamilyIndices Renderer::findQueueFamilies(VkPhysicalDevice device)
{
	QueueFamilyIndices indices;

	std::uint32_t queueFamilyCount = 0;
	vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, nullptr);

	std::vector<VkQueueFamilyProperties> queueFamilies(queueFamilyCount);
	vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, queueFamilies.data());

	std::uint32_t i = 0;
	for (const auto& queueFamily : queueFamilies)
	{
		if (queueFamily.queueFlags & VK_QUEUE_GRAPHICS_BIT)
		{
			indices.graphics = i;
		}

		VkBool32 presentSupport = false;
		vkGetPhysicalDeviceSurfaceSupportKHR(device, i, _vkSurface, &presentSupport);
		if (presentSupport)
		{
			indices.present = i;
		}

		// Stop the search if all indices required are found
		if (indices.isComplete())
		{
			break;
		}
		++i;
	}

	return indices;
}

void Renderer::createLogicalDevice()
{
	QueueFamilyIndices indices = findQueueFamilies(_vkPhysicalDevice);

	std::vector<VkDeviceQueueCreateInfo> queueCreateInfos;
	std::set<std::uint32_t> uniqueQueueFamilies =
	{
		indices.graphics.value(),
		indices.present.value()
	};
	float queuePriority = 1.0f;
	for (const auto& queueFamily : uniqueQueueFamilies)
	{
		VkDeviceQueueCreateInfo queueCreateInfo{};
		queueCreateInfo.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
		queueCreateInfo.pNext = nullptr;
		queueCreateInfo.flags = 0;
		queueCreateInfo.queueFamilyIndex = indices.graphics.value();
		queueCreateInfo.queueCount = 1;
		queueCreateInfo.pQueuePriorities = &queuePriority;
		queueCreateInfos.emplace_back(queueCreateInfo);
	}

	VkPhysicalDeviceFeatures deviceFeatures{};

	VkDeviceCreateInfo createInfo{};
	createInfo.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
	createInfo.pNext = nullptr;
	createInfo.flags = 0;
	createInfo.queueCreateInfoCount = static_cast<std::uint32_t>(queueCreateInfos.size());
	createInfo.pQueueCreateInfos = queueCreateInfos.data();
	createInfo.enabledExtensionCount = 0;
	createInfo.ppEnabledExtensionNames = nullptr;
	createInfo.pEnabledFeatures = &deviceFeatures;
	if (enableValidationLayers)
	{
		createInfo.enabledLayerCount = static_cast<std::uint32_t>(vkValidationLayers.size());
		createInfo.ppEnabledLayerNames = vkValidationLayers.data();
	}
	else
	{
		createInfo.enabledLayerCount = 0;
	}

	if (vkCreateDevice(_vkPhysicalDevice, &createInfo, nullptr, &_vkDevice) != VK_SUCCESS)
	{
		ERROR("Failed to create logical device.")
	}

	vkGetDeviceQueue(_vkDevice, indices.graphics.value(), 0, &_vkGraphicsQueue);
	vkGetDeviceQueue(_vkDevice, indices.present.value(), 0, &_vkPresentQueue);
}

void Renderer::createSurface()
{
	_window->windowCreateSurface(_vkInstance, &_vkSurface);
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
	createInfo.flags = 0;
	createInfo.pApplicationInfo = &appInfo;
	VkDebugUtilsMessengerCreateInfoEXT debugCreateInfo;
	if (enableValidationLayers)
	{
		createInfo.enabledLayerCount = static_cast<std::uint32_t>(vkValidationLayers.size());
		createInfo.ppEnabledLayerNames = vkValidationLayers.data();
		populateDebugMessengerCreateInfo(debugCreateInfo);
		createInfo.pNext = (VkDebugUtilsMessengerCreateInfoEXT*)&debugCreateInfo;
	}
	else
	{
		createInfo.enabledLayerCount = 0;
		createInfo.ppEnabledLayerNames = nullptr;
		createInfo.pNext = nullptr;
	}
	std::vector<const char*> vkExtensions = getVkRequiredExtensions();
	createInfo.enabledExtensionCount = static_cast<std::uint32_t>(vkExtensions.size());
	createInfo.ppEnabledExtensionNames = vkExtensions.data();

	if (vkCreateInstance(&createInfo, nullptr, &_vkInstance))
	{
		ERROR("Cannot create Vulkan instance.");
	}
}

std::vector<const char*> Renderer::getVkRequiredExtensions()
{
	std::pair<const char**, std::uint32_t> windowRequiredInstanceExtensions = _window->windowGetRequiredInstanceExtensions();
	std::vector<const char*> vkExtensions(windowRequiredInstanceExtensions.first, windowRequiredInstanceExtensions.first + windowRequiredInstanceExtensions.second);
	if (enableValidationLayers)
	{
		vkExtensions.emplace_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
	}
	return vkExtensions;
}

void Renderer::setupDebugMessenger()
{
	if (!enableValidationLayers)
	{
		return;
	}
	VkDebugUtilsMessengerCreateInfoEXT createInfo;
	populateDebugMessengerCreateInfo(createInfo);

	PFN_vkCreateDebugUtilsMessengerEXT vkCreateDebugUtilsMessenger = (PFN_vkCreateDebugUtilsMessengerEXT)vkGetInstanceProcAddr(_vkInstance, "vkCreateDebugUtilsMessengerEXT");
	if (vkCreateDebugUtilsMessenger == nullptr)
	{
		ERROR("Unable to get vkCreateDebugUtilsMessengerEXT function pointer.");
	}
	vkCreateDebugUtilsMessenger(_vkInstance, &createInfo, nullptr, &_debugMessenger);
}

void Renderer::populateDebugMessengerCreateInfo(VkDebugUtilsMessengerCreateInfoEXT& createInfo)
{
	createInfo = {};
	createInfo.sType = VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT;
	createInfo.pNext = nullptr;
	createInfo.flags = 0;
	createInfo.messageSeverity = VK_DEBUG_UTILS_MESSAGE_SEVERITY_VERBOSE_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT;
	createInfo.messageType = VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT;
	createInfo.pfnUserCallback = vkDebugCallback;
	createInfo.pUserData = nullptr;
}

void Renderer::destroyDebugMessenger()
{
	PFN_vkDestroyDebugUtilsMessengerEXT vkDestroyDebugUtilsMessenger = (PFN_vkDestroyDebugUtilsMessengerEXT)vkGetInstanceProcAddr(_vkInstance, "vkDestroyDebugUtilsMessengerEXT");
	if (vkDestroyDebugUtilsMessenger == nullptr)
	{
		ERROR("Unable to get vkDestroyDebugUtilsMessengerEXT function pointer.");
	}
	vkDestroyDebugUtilsMessenger(_vkInstance, _debugMessenger, nullptr);
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
