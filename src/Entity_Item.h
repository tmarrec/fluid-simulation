#pragma once

#include <memory>

#include "Entity.h"

class Entity_Item {

public:
	Entity_Item(std::unique_ptr<Entity> & shape_ptr);
	~Entity_Item(void);
	void prout();

	const std::string name() const;

private:
	const std::unique_ptr<Entity> & _shape_ptr;
	const std::string _name;
};
