#include "ECS.h"

#include "CameraComponent.h"

Entity::Entity(MsgBus_ptr messageBus)
	: System{messageBus}
{
	_id = _counter++;
	Message helloMsg {HELLO, this};
	postMessage(helloMsg);
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

void Entity::cout(std::string string) const
{
	std::cout << "0x" << std::hex << std::this_thread::get_id() << " ";
	std::cout << "  \033[100m\033[1m";
	std::cout << "[Entity " << getEntityID() << "]";
	std::cout << "\033[49m\033[0m";
	std::cout << " " << string << std::endl;
}

void Entity::entityPostMessage(Message & msg)
{
	postMessage(msg);
}
		
void Entity::handleMessage(Message & msg)
{
	switch(msg._type)
	{
		case HELLO_ACK:
			cout("Loaded in the \033[45m\033[1m[MessageBus]\033[49m\033[0m");
			break;
			
		case ASK_CAMERA_INFOS_FOR_DRAW:
			if (hasComponent<CameraComponent>())
			{
				auto cameraComponent = getComponent<CameraComponent>();
				auto view = cameraComponent.view();
			}
			break;

		case GL_SIZE:
			if (hasComponent<CameraComponent>())
			{
				auto cameraComponent = getComponent<CameraComponent>();
				cameraComponent.setProjection(msg._width, msg._height);
			}
			break;
	
		default:
			break;
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


ID Entity::_counter = 0;
ID Component::_counter = 0;

Manager::Manager(MsgBus_ptr messageBus)
	: System{messageBus}
{
	Message helloMsg {HELLO, this};
	postMessage(helloMsg);
}

void Manager::cout(std::string string) const
{
	std::cout << "0x" << std::hex << std::this_thread::get_id() << " ";
	std::cout << "       \033[42m\033[1m";
	std::cout << "[ECS]";
	std::cout << "\033[49m\033[0m";
	std::cout << " " << string << std::endl;
}
		
void Manager::handleMessage(Message & msg)
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

Entity & Manager::addEntity(MsgBus_ptr messageBus)
{
	Entity * e = new Entity(messageBus);
	std::unique_ptr<Entity> uPtr {e};
	_entities.emplace_back(std::move(uPtr));
	return *e;
}

