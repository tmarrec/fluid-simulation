#include "Entity_Item.h"
#include <iostream>

Entity_Item::Entity_Item(std::shared_ptr<Entity> shape_ptr)
	: _shape_ptr{shape_ptr}
	, _name{shape_ptr->name()}
{
	
}

Entity_Item::~Entity_Item(void) {

}

const std::string Entity_Item::name() const {
	return _name;
}

std::shared_ptr<Entity> Entity_Item::shape_ptr() {
	return _shape_ptr;
}

