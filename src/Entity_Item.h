#pragma once

#include <memory>

#include "Entity.h"

class Entity_Item {

public:
	Entity_Item(std::shared_ptr<Entity> entity_ptr);
	~Entity_Item(void);

	const std::string name() const;
	std::shared_ptr<Entity> entity_ptr();

private:
	const std::shared_ptr<Entity> _entity_ptr;
	const std::string _name;
};
