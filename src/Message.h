#pragma once

class System;

enum Type
{
	 TEST
	,HELLO
	,HELLO_ACK
	,INIT_GL
	,DRAW
	,RESIZE_GL
};

class Message
{
public:
	Message(Type type);
	Message(Type type, System * system);
	Message(Type type, int width, int height);
	const Type & type() const;
	System * system() const;
	int width() const;
	int height() const;

private:
	const Type _type;
	System * _system = nullptr;
	int _width = 0;
	int _height = 0;
};
