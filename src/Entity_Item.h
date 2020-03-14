#pragma once

#include <memory>

#include "Entity.h"

class Entity_Item {

public:
	Entity_Item(std::shared_ptr<Entity> shape_ptr);
	~Entity_Item(void);

	const std::string name() const;
	std::shared_ptr<Entity> shape_ptr();

private:
	const std::shared_ptr<Entity> _shape_ptr;
	const std::string _name;
};