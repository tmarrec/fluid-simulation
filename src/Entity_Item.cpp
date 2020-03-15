#include "Entity_Item.h"
#include <iostream>

Entity_Item::Entity_Item(std::shared_ptr<Entity> entity_ptr)
	: _entity_ptr{entity_ptr}
{
	switch(entity_ptr->type()) {
		case CAMERA:
			_name = "Camera";
			break;
		case LIGHT:
			_name = "Light";
			break;
		case SHAPE:
			_name = "Shape";
			break;
	}
	_name += " "+std::to_string(entity_ptr->id());
}

Entity_Item::~Entity_Item(void) {

}

std::string Entity_Item::name() const {
	return _name;
}

std::shared_ptr<Entity> Entity_Item::entity_ptr() {
	return _entity_ptr;
}

