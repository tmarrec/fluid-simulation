#pragma once

class System;

enum Type
{
	TEST
	,HELLO
	,HELLO_ACK
};

class Message
{
public:
	Message(Type type);
	Message(Type type, System * system);
	const Type & type() const;
	const System * system() const;

private:
	const Type _type;
	System * _system;
};
