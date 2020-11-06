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
	/* temp circle motion */
	float angle = 0;
	float speed = (2*M_PI)/5;
	float radius = 100;
	void update([[maybe_unused]] double _deltaTime) override
	{
		radius = _intensity*100;
		speed = (2*M_PI)/_intensity*0.2;
		angle += speed*_deltaTime;
		auto pos = entity->getComponent<TransformComponent>().position();
		pos.x = cos(angle)*radius;
		pos.z = sin(angle)*radius;
		entity->getComponent<TransformComponent>().setPosition(pos);
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
	glm::mat4 getLightSpaceMatrix() const
	{
		// TODO add near plane and far plane ------------------------------->~~~~~~~~~~~~~~~~
		glm::mat4 lightProjection = glm::ortho(-10.0f, 10.0f, -10.0f, 10.0f, 0.0f, 7000000.5f);
		glm::mat4 lightView = glm::lookAt(getPosition(), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
		return lightProjection*lightView;
	}


private:
	Renderer__ _renderer;
	glm::vec3 _color;
	double _intensity;
};
