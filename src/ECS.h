#pragma once

#include <bitset>
#include <array>

class Component;
class Entity;

using ComponentID = std::uint64_t;

inline ComponentID getComponentTypeID()
{
	static ComponentID lastID = 0;	
	return lastID++;
}

template <typename T> inline ComponentID getComponentTypeID() noexcept
{
	static ComponentID typeID = getComponentTypeID();
	return typeID;
}

constexpr std::uint64_t maxComponents = 64;

using ComponentBitSet = std::bitset<maxComponents>;
using ComponentArray = std::array<Component *, maxComponents>;

class Component
{
public:
	Entity* entity;

	Component()
	{
		_id = _counter++;
	}

	ComponentID getComponentID() const
	{
		return _id;
	}

	virtual void init() = 0;
	virtual void update() = 0;
	virtual void draw()
	{

	}

	virtual ~Component()
	{

	}
	static ComponentID _counter;

private:
	ComponentID _id;

};

ComponentID Component::_counter = 0;

class Entity
{
public:
	void update()
	{
		//std::vector<Message> messages;
		for(const auto & component : _components)
		{
			component->update();
		}
		for(const auto & component : _components)
		{
			//auto msg = component->draw();
			//messages.emplace_back(msg);
			component->draw();
		}
	}

	void draw()
	{

	}

	bool isActive() const
	{
		return _active;
	}

	void destroy()
	{
		_active = false;
	}

	template <typename T> bool hasComponent() const
	{
		return _componentBitSet[getComponentTypeID<T>];
	}

	template <typename T, typename... TArgs>
	T& addComponent(TArgs&&... mArgs)
	{
		T* c {new T(std::forward<TArgs>(mArgs)...)};
		c->entity = this;
		std::unique_ptr<Component> uPtr {c};
		_components.emplace_back(std::move(uPtr));

		auto componentTypeID = getComponentTypeID<T>();
		_componentArray[componentTypeID] = c;
		_componentBitSet[componentTypeID] = true;

		c->init();
		return *c;
	}

	template <typename T> T& getComponent() const
	{
		auto ptr {_componentArray[getComponentTypeID<T>()]};
		return *static_cast<T*>(ptr);
	}


private:
	bool _active = true;
	std::vector<std::unique_ptr<Component>> _components;

	ComponentArray _componentArray;
	ComponentBitSet _componentBitSet;

};

class Manager : public System
{
public:
	Manager(MsgBus_ptr messageBus)
	: System{messageBus}
	{
		Message helloMsg {HELLO, this};
		postMessage(helloMsg);
	}

	void cout(std::string string) const
	{
		std::cout << "0x" << std::hex << std::this_thread::get_id() << " ";
		std::cout << "       \033[42m\033[1m";
		std::cout << "[ECS]";
		std::cout << "\033[49m\033[0m";
		std::cout << " " << string << std::endl;
	}
		
	void handleMessage(Message & msg)
	{
		switch(msg._type)
		{
			case HELLO_ACK:
				cout("Loaded in the \033[45m\033[1m[MessageBus]\033[49m\033[0m");
				break;
	
			default:
				break;
		}
	}

	void update()
	{
		for (const auto & entity : _entities)
		{
			entity->update();
		}
	}

	void draw()
	{
		for (const auto & entity: _entities)
		{
			entity->draw();
		}
	}

	void refresh()
	{
		_entities.erase(std::remove_if(std::begin(_entities), std::end(_entities),
			[](const std::unique_ptr<Entity> &mEntity)
			{
				return !mEntity->isActive();
			}),
			std::end(_entities));
	}

	Entity & addEntity()
	{
		Entity * e = new Entity();
		std::unique_ptr<Entity> uPtr {e};
		_entities.emplace_back(std::move(uPtr));
		return *e;
	}

private:
	std::vector<std::unique_ptr<Entity>> _entities;

};
