#pragma once


#include "ECS.h"
#include "Renderer.h"
#include "TransformComponent.h"
#include "glm/glm/ext/matrix_clip_space.hpp"
#include "glm/glm/ext/matrix_transform.hpp"
#include <cstdint>

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
		radius = _intensity*30;
		speed = (2*M_PI)/_intensity*0.005;
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
	std::vector<glm::mat4> getLightSpaceMatrices(std::uint64_t __depthShadowWidth, std::uint64_t __depthShadowHeight) const
	{
		glm::mat4 lightProjection = glm::infinitePerspective(glm::radians(90.0f), (float)__depthShadowWidth/(float)__depthShadowHeight, 0.1f);
		std::vector<glm::mat4> lightSpaceMatrices;
		auto lightPos = getPosition();
		lightSpaceMatrices.emplace_back(lightProjection*glm::lookAt(lightPos, lightPos+glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f)));
		lightSpaceMatrices.emplace_back(lightProjection*glm::lookAt(lightPos, lightPos+glm::vec3(-1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f)));
		lightSpaceMatrices.emplace_back(lightProjection*glm::lookAt(lightPos, lightPos+glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f)));
		lightSpaceMatrices.emplace_back(lightProjection*glm::lookAt(lightPos, lightPos+glm::vec3(0.0f, -1.0f, 0.0f), glm::vec3(0.0f, 0.0f, -1.0f)));
		lightSpaceMatrices.emplace_back(lightProjection*glm::lookAt(lightPos, lightPos+glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f, -1.0f, 0.0f)));
		lightSpaceMatrices.emplace_back(lightProjection*glm::lookAt(lightPos, lightPos+glm::vec3(0.0f, 0.0f, -1.0f), glm::vec3(0.0f, -1.0f, 0.0f)));
		return lightSpaceMatrices;
	}


private:
	Renderer__ _renderer;
	glm::vec3 _color;
	double _intensity;
};
