#include "Entity_Item.h"
#include <iostream>

Entity_Item::Entity_Item(std::unique_ptr<Entity> & shape_ptr)
	: _shape_ptr{shape_ptr}
	, _name{std::to_string(shape_ptr->id())}
{
	
}

Entity_Item::~Entity_Item(void) {

}

const std::string Entity_Item::name() const {
	return _name;
}
