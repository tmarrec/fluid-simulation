#include "ECS.h"

#include "CameraComponent.h"
#include "src/DrawableComponent.h"
#include <algorithm>

Entity::Entity()
{
	_id = _counter++;
}

ID Entity::getEntityID() const
{
	return _id;
}

void Entity::update()
{
	for(const auto & component : _components)
	{
		component->update();
	}
	for(const auto & component : _components)
	{
		component->draw();
	}
}

void Entity::init()
{
	std::cout << "teste" << std::endl;
}

void Entity::draw()
{

}

bool Entity::isActive() const
{
	return _active;
}

void Entity::destroy()
{
	_active = false;
}


ID Entity::_counter = 0;
ID Component::_counter = 0;

void Manager::update()
{
	for (const auto & entity : _entities)
	{
		entity->update();
	}
}

void Manager::draw()
{
	for (const auto & entity: _entities)
	{
		entity->draw();
	}
}

void Manager::refresh()
{
	_entities.erase(std::remove_if(std::begin(_entities), std::end(_entities),
		[](const std::unique_ptr<Entity> &mEntity)
		{
			return !mEntity->isActive();
		}),
		std::end(_entities));
}

Entity & Manager::addEntity()
{
	Entity * e = new Entity();
	std::unique_ptr<Entity> uPtr {e};
	_entities.emplace_back(std::move(uPtr));
	return *e;
}

