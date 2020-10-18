#pragma once


#include "ECS.h"
#include "Renderer.h"
#include "TransformComponent.h"

using Renderer__ = std::shared_ptr<Renderer>; 

class LightComponent : public Component
{
public:
	LightComponent(Renderer__ __renderer, glm::vec3 __color, double __intensity)
	: Component{}
	, _renderer { __renderer }
	, _color { __color }
	, _intensity { __intensity }
	{}

	void init() override
	{
		ASSERT(entity->hasComponent<TransformComponent>(), "entity should have a TransformComponent");
		_renderer->addLight(this);
	}
	void draw() override
	{
	}
	void update([[maybe_unused]] double _deltaTime) override
	{
		entity->getComponent<TransformComponent>().move({0.02f, 0.0f, 0.02f});
	}
	~LightComponent() override
	{
	}

	glm::vec3 color() const { return _color; }
	double intensity() const { return _intensity; }
	glm::vec3 getPosition() const
	{
		ASSERT(entity->hasComponent<TransformComponent>(), "entity should have a TransformComponent");
		return entity->getComponent<TransformComponent>().position();
	}


private:
	Renderer__ _renderer;
	glm::vec3 _color;
	double _intensity;
};
