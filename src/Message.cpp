#include "Message.h"

#include "System.h"

Message::Message(Type type)
	: _type{type}
{

}

Message::Message(Type type, System * system)
	: _type{type}
	, _system{system}
{

}

Message::Message(Type type, int width, int height)
	: _type{type}
	, _width{width}
	, _height{height}
{

}

const Type & Message::type() const
{
	return _type;
}

System * Message::system() const
{
	return _system;
}

int Message::width() const
{
	return _width;
}

int Message::height() const
{
	return _height;
}
