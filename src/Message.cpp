#include "Message.h"

#include "System.h"

Message::Message(Type type)
	: _type{type}
	, _system{nullptr}
{

}

Message::Message(Type type, System * system)
	: _type{type}
	, _system{system}
{

}

const Type & Message::type() const
{
	return _type;
}

const System * Message::system() const
{
	return _system;
}
