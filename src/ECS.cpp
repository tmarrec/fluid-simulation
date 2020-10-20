#include "ECS.h"

#include "CameraComponent.h"
#include "src/DrawableComponent.h"
#include <algorithm>
#include <execution>
#include <pstl/glue_execution_defs.h>

ID Entity::_counter = 0;
ID Component::_counter = 0;

Entity::Entity()
{
	_id = _counter++;
}

ID Entity::getEntityID() const
{
	return _id;
}

void Entity::update(double __deltaTime)
{
	for(const auto & component : _components)
	{
		component->update(__deltaTime);
	}
}

void Entity::init()
{

}

void Entity::draw()
{
	for(const auto & component : _components)
	{
		component->draw();
	}
}

bool Entity::isActive() const
{
	return _active;
}

void Entity::destroy()
{
	_active = false;
}

void ECS_Manager::update(double __deltaTime)
{
	refresh();
	for (const auto & entity: _entities)
	{
		entity->update(__deltaTime);	
	}
}

void ECS_Manager::draw()
{
	for (const auto & entity: _entities)
	{
		entity->draw();
	}
}

void ECS_Manager::refresh()
{
	_entities.erase(std::remove_if(std::begin(_entities), std::end(_entities),
		[](const std::unique_ptr<Entity> &mEntity)
		{
			return !mEntity->isActive();
		}),
		std::end(_entities));
}

Entity & ECS_Manager::addEntity()
{
	Entity * e = new Entity();
	std::unique_ptr<Entity> uPtr {e};
	uPtr->init();
	_entities.emplace_back(std::move(uPtr));
	return *e;
}

