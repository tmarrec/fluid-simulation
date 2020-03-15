#include "Entity_Item.h"
#include <iostream>

Entity_Item::Entity_Item(std::shared_ptr<Entity> entity_ptr)
	: _entity_ptr{entity_ptr}
	, _name{entity_ptr->name()}
{
	
}

Entity_Item::~Entity_Item(void) {

}

const std::string Entity_Item::name() const {
	return _name;
}

std::shared_ptr<Entity> Entity_Item::entity_ptr() {
	return _entity_ptr;
}

